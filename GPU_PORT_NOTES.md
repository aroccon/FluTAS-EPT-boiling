# GPU port on Leonardo: bugs found and fixed (19–20 Oct 2026)

Case: `src/single_nuc_site` (32×80×32, IBM wall with `IBM_h = 4`, y-outflow, SEM inflow,
Smagorinsky SGS, phase change), app `two_phase_cdi_marangoni`, restart from step 1680000.

**Result:** the 4-GPU run (Leonardo Booster, A100, NVHPC 25.3) now reproduces the CPU run
(GNU build, 4 MPI ranks). Field checksums agree to 6–7 digits after every stage of the time
step. Before these fixes the GPU run blew up to NaN within 5 time steps.

## How the bugs were found

1. NaN/Inf counts per field after each stage of the time step (`chk_fields`).
2. Largest absolute value per field, to catch blow-ups before they become NaN.
3. Interior checksums (sum of |value|) per field, compared line by line between CPU and GPU.
4. A y-halo consistency check against the neighbour ranks (plain MPI on the host).
5. Single-GPU vs 4-GPU runs, to separate GPU-kernel bugs from decomposition bugs.

## Bugs and fixes

### 1. GPU halo exchange overwrote physical boundaries (`bound.f90`, `updthalo`)
- **Symptom:** NaN from step 1; `w(i,80,0) = -huge`, which became Inf through the y-outflow BC.
- **Cause:** on GPU, `updthalo` fills its receive buffers with `-huge` and copies them into the
  halo unconditionally. With `MPI_PROC_NULL` (a physical boundary) nothing is received, so
  the BC values that `bounduvw` had set just before were replaced with `-huge`. The CPU path
  (`MPI_SENDRECV` on derived types) leaves the halo untouched, which is why CPU runs were fine.
- **Fix:** unpack a received buffer only if that neighbour exists (x, y and z).

### 2. IBM wall requires z not to be split (`targets/target.generic-gpu`)
- **Symptom:** with z split over 4 GPUs, `w` differed from the CPU after the momentum step.
- **Cause:** about 190 uses of `IBM_h` in 15 files treat the local `k <= IBM_h` as the solid,
  so every rank put a wall in its own bottom 4 layers.
- **Fix:** GPU build uses `-D_DECOMP_Z` (as the CPU build already did): z is never split,
  only x and y via `dims_in`. With `dims_in = 1 4`, y is split in 4.
- **Limitation:** `nx`, `ny` must be divisible by `dims_in`; z cannot be split.

### 3. `dl` was a static managed array, so it was 0 on GPUs 1–3 (`param.f90`, `di_acdi.f90`, main)
- **Symptom:** `mu_t = 0` on ranks 1–3 only; rank 0 correct.
- **Cause:** `dl(3)` was fixed-size with `!@cuf attributes(managed)`. Static managed variables
  are created at program start on GPU 0, before each rank selects its GPU, so kernels on
  GPUs 1–3 read zeros. (`initmpi.f90` already warns that all GPU memory must be dynamic.)
- **Fix:** removed the managed attribute from `dl` in `param.f90` and from the `dl` argument
  of `rk_ent_ibm` in `di_acdi.f90`; added `!$acc enter data copyin(dl)` in main right after
  `initmpi`. `dl` was the only static managed module variable in the code.

### 4. Poisson eigenvalues prepared for the wrong solver (`initsolver.f90`)
- **Symptom:** pressure about 10% off on 1 GPU; wrong velocity after projection on 4 GPUs.
- **Cause:** under `_OPENACC`, `initsolver` builds `lambdaxy` for `solver_gpu`: cuFFT R2C
  ordering for periodic directions, and a transposed z-pencil shape and rank indexing. But
  this app always calls `solver_cpu` (FFTW r2r + IBM), also in the GPU build.
- **Fix:** the GPU layout is used only with `-D_USE_SOLVER_GPU` (not set). Otherwise the CPU
  layout is used.

### 5. Smaller GPU fixes
- `sgs.f90`:
  - GPU loops used the host array `n(:)` with `default(present)` → scalar copies `n1,n2,n3`.
  - `is_bound` was set only on the first call but not `save`d → now `save`.
  - SGS kernels run on `async(1)` but nothing waited for them → `!$acc wait(1)` at the end.
- `fillps.f90`: `deallocate(outflow_ind)` was inside a `!$acc kernels` region → moved after it.
- `chkdt.f90` (`chkdt_tw`): `max_visc`, `max_cond` added to the GPU `reduction` clause.
- `targets/target.generic-gpu`:
  - `-gpu=cc80` (A100) instead of `cc70`, via `GPU_CC ?= 80`.
  - `mem:managed`, so the ~40 allocatable arrays without a managed attribute (`mu_t`,
    `rhocp`, `psiold`, …) are also on the GPU.
  - FFTW rpath.
  - `GPU_MPI ?= 1` switch (`GPU_MPI=0` keeps MPI buffers in managed memory). The MPI
    library is CUDA-aware, and both settings give identical results.
- `initmpi.f90`: `acc_set_device_num` after `cudaSetDevice`, and the startup line prints both
  devices. Not the cause of any bug here, but recommended for MPI + OpenACC.

## Input and build fixes

- `sem.f90`: `neddy = 2000` to match `single_nuc_site/sem_restart_1.txt`, which holds 2000 eddies.
  The code had 150000, so the restart aborted.
- `single_nuc_site/dns.in`: `dims_in = 1 4`. These files are git-ignored, so copy them by hand.
- `make.deps`: 13 missing module dependencies added, so `make -j` and incremental builds are safe.
- `Makefile`: `ln -sf` so incremental builds don't fail on the existing `flutas` link.
- `compile_leo_gpu.sh`: `make -j $NJOBS`.
- `compile_cpu.sh`: `module purge`, Leonardo modules, Leonardo FFTW path.
- `single_nuc_site/go_leo.sh`: Leonardo GPU job (1 rank per GPU, `binder_leo.sh`). The `fftw`
  module was removed because it brought in a second MPI.
- `single_nuc_site/go_leo_cpu.sh`: CPU reference job.
- `binder_leo.sh`: one network card per rank.

## Debug tools

- `debug.f90`: `chk_fields(label, u,v,w,p,psi,ent,mu_t,rho,mu)` and `chk_arr(label, a)` print
  NaN/Inf counts (`NANCHECK`), maxima (`NANMAX`), interior checksums (`NANSUM`/`ARRSUM`) and
  y-halo mismatches (`YHALO`). They do nothing unless the environment variable
  `FLUTAS_NANCHECK` is set.
- The calls were removed after the fixes (21 Oct). To debug again, add
  `call chk_fields('label',u,v,w,p,psi,ent,mu_t,rho,mu)` after the stages of interest in the main
  program, and uncomment `export FLUTAS_NANCHECK=1` in `go_leo.sh` / `go_leo_cpu.sh`.

## Open items

1. Run the 4-GPU case for longer and check that `dt` grows as on the CPU (6.4e-7 limit).
2. `sgs.f90` (van Driest): `dw(6) = dl(3)*(n3-j+0.5)` should probably use `k`. It doesn't
   matter here, because only the bottom z wall is a wall.
3. `neddy` is hard-coded. Branch `lw_clean_examples` reads it from `sem_input_*.txt`. Worth
   merging, together with its tagging guard (`do_tagging`) and its `fillps` changes.
   Do **not** take its `cmpt_umax` debug prints or its hard-coded `u_reg = 0.1`.
4. `load.f90` doesn't check whether `MPI_FILE_OPEN` succeeded: a missing checkpoint crashes in
   `MPI_File_get_size` instead of printing an error.
5. Performance: arrays that are never prefetched (`mu_t`, `rhocp`, …) and per-step managed
   allocations (`outflow_ind` in `fillps`) may cost time. Check with `nsys`.
6. Making the IBM code decomposition-aware (offset `IBM_h` by `ijk_start(3)`) would allow splitting z.
