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
  - `GPU_MPI ?= 1` switch. Note: `GPU_MPI=0` does **not** compile with `_EPA2A`, so only
    `GPU_MPI=1` is usable. The MPI library (HPC-X 2.19) is CUDA-aware.
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

## Performance (21 Oct, 300 steps from step 1680000, profiler report)

| section | CPU, 4 cores (GNU) | GPU, 4×A100 |
|---|---|---|
| total per step | 13 ms | 836 ms (~63× slower) |
| `SOLVER` (300 steps) | 0.89 s | 237 s |
| 6 transposes (300 steps) | 0.18 s (~0.03 ms/call) | 228 s (25–65 ms/call) |
| everything else (300 steps) | 3.1 s | 14 s |

- The GPU build always calls `solver_cpu` (FFTW + host transposes + `gaussel`), because only it
  handles the IBM wall. `solver_gpu` is never called.
- The transposes are ~1300× slower in the GPU build than on CPU, with the same code and data size.
  Likely cause: the CUDA-aware MPI stack (HPC-X/UCX) handling `MPI_Alltoallw` with derived types on
  host or managed memory. Tried (21 Oct): `px`, `py`, `pz` in `solver_cpu` as pinned host arrays,
  allocated once. Transposes stayed at ~29 ms and the run went to NaN at step 4, so it was
  reverted. Next: `nsys` profile of a few steps to see what the MPI calls actually do.
- **Fix 1 (21 Oct):** every 2decomp transpose did `decomp = decomp_main`, a copy of a derived
  type with allocatable components. Under `mem:managed` each component copy is a synchronous
  `cudaMemcpy` (seen with `nsys --cudabacktrace`: `transpose_*_real` → `pgf90_man_copyover`).
  Now `decomp` is a pointer (`decomp_main` and `opt_decomp` are `target`). STEP for 300 steps:
  250 s → 98 s. x↔y transposes 26 → 0.8 s, y↔z 29 → 3.2 s; z↔x still 34 s.
- **Fix 2 (21 Oct):** in the GPU build, `solver_cpu` does z↔x via y (two `Alltoallv`
  transposes) instead of the direct `MPI_Alltoallw` with derived types, which the CUDA-aware MPI
  handles very slowly on managed memory. The `contiguous` attribute added to the host transposes
  earlier had no measurable effect.
- **Result after both fixes** (300 steps, 4×A100, clean run): STEP 250 s → **40 s** (0.13 s/step,
  6× faster); SOLVER 237 s → 26.5 s, of which 17.3 s in transposes (y↔z still ~3.5 ms/call vs
  0.03 ms on CPU). The rest of the step is ~13.5 s (SEM 3.3, RK 3.1, PC 2.3, HT 1.9, DI 1.1).
  The CPU is still ~10× faster on this tiny grid (13 ms/step).
- The rest is ~4.5× slower than CPU: managed-memory pages move between CPU loops (`bound.f90`, most
  of `di_acdi.f90`) and GPU kernels at every step.
- This case is far too small for GPUs (82k cells, ~20k per GPU). Measure performance on a
  realistic grid.
- Every run overwrites `data/restart_dir` **and** `sem_restart_1.txt` at the end. Keep a backup of
  the starting state (`checkpoints/single_nuc_site_1680000` on Leonardo).

## GPU pressure solver (22 Oct)

- Switch: `USE_SOLVER_GPU` (default 1 = `solver_gpu`; `USE_SOLVER_GPU=0 ./compile_leo_gpu.sh` = `solver_cpu`).
- `fft.f90`: `ND` (cell-centred) pressure BC on GPU as DCT-IV (FFTW `REDFT11`), built on the existing
  DCT-II/III: forward = scale by 2cos(pi(2n+1)/(4N)), DCT-II, recurrence Y_0 = W_0/2, Y_k = W_k - Y_(k-1);
  backward = transposed chain. Checked against direct transforms (errors ~N·eps).
- `solver_gpu.f90`: IBM offset in the tridiagonal solve (`k = 1+IBM_h ... nz`), as in `solver_cpu`.
- Validation: identical to `solver_cpu` on 1 GPU and on 4 GPUs (right-hand side and solution of all
  3 solves per step, 3 steps, 7 digits).
- Also fixed: `solver_gpu` called `cudaMemAdvise` on the device arrays `pxc`, `pyc_t`, `py_t`. That
  fails with 'invalid argument' and left the error behind, so UCC's `MPI_Allreduce` on GPU buffers
  then reported `cudaGetLastError() failed` and aborted the reduction (4 GPUs only).
- Performance, 300 steps, 4×A100: STEP 40 s → **20 s** (0.067 s/step); SOLVER 26.5 s → 6.7 s;
  transposes 17.3 s → 0.54 s. Overall 0.84 → 0.067 s/step since 21 Oct (12×). The rest of the
  step (~13 s: SEM 3.3, RK 3.1, PC 2.3, HT 1.7, DI 1.2) is now the larger part.
- 2000 steps from step 1680000, 4 GPUs, `solver_gpu` vs `solver_cpu`: `Max T_wall` and `dt` agree to
  ~7 digits at the end (differences grow slowly from ~1e-9, normal rounding growth). STEP 133 s vs
  268 s. **`USE_SOLVER_GPU=1` is now the default** (target and `compile_leo_gpu.sh`).
- `restore_case.sh <case> <nstep>` restores the clean starting state before comparisons.

## Porting CPU routines to the GPU (23 Oct)

Profiling (`nsys` with page-fault tracing) showed that inside a step the CPU is busy all the time
and the GPU is mostly idle: many CPU loops touch managed arrays, so data moves back and forth
(~32 MB each way per step per rank). Porting single routines gains little (the faults move to the
next CPU routine), so routines are ported in batches, each validated over 2000 steps against the
previous version (`Max T_wall`, `dt`: 13–14 digits early, ~9–11 digits after 2000 steps).

| batch | routines | STEP (2000 steps) |
|---|---|---|
| start | (`solver_gpu`) | 133.3 s |
| SEM | `SEM_iso_dir` eddy loop (one thread per inlet point) | 133.5 s |
| 1 | `make_vol_avg_vel`, `mom_to_vel`, `set_IBM_vel_to_zero`, `add_a_to_b`, `cmpt_u_grad_rho`, `update_property_IBM`, `mfx/mfy/mfz` loops in main | 129.4 s |
| 2 | `pres_tw_src`, `copy_fields`, `psi_to_ls`, `cmpt_umax` (max reduction), `cmpt_norm`, `SEM_apply` | **106.8 s** (53 ms/step) |

Rules used: `!$acc parallel loop collapse(3) default(present)` on managed arrays; inside kernels only
local scalars for `n(:)`, `IBM_h`, `rho1-rho2`, `dl`, `dli` (module arrays like `dli` are static, not
on the GPU). Profiler sub-regions `RK:*` and `SEM:*` added (profiler limit raised to 48 labels).
Remaining per step: SOLVER 11 ms, SEM 9.4 ms, RK 8.9 ms, HT 7.3 ms, PC 5.6 ms, DI 5.2 ms.
Comparison with the CPU build (GNU, 4 MPI ranks on 4 cores of one node, same case): 13 ms/step.
The 4-GPU run (53 ms/step) is still ~4× slower on this tiny grid (~20k cells per GPU); a full CPU
node (32 cores) would be faster still. GPU gains have to be measured on a production-size grid.

## Open items

1. Run the 4-GPU case for longer and check that `dt` grows as on the CPU (6.4e-7 limit).
2. `sgs.f90` (van Driest): `dw(6) = dl(3)*(n3-j+0.5)` should probably use `k`. It doesn't
   matter here, because only the bottom z wall is a wall.
3. `neddy` is hard-coded. Branch `lw_clean_examples` reads it from `sem_input_*.txt`. Worth
   merging, together with its tagging guard (`do_tagging`) and its `fillps` changes.
   Do **not** take its `cmpt_umax` debug prints or its hard-coded `u_reg = 0.1`.
4. `load.f90` doesn't check whether `MPI_FILE_OPEN` succeeded: a missing checkpoint crashes in
   `MPI_File_get_size` instead of printing an error.
5. Performance: fix the slow transposes first (see above). Then move the remaining CPU loops to
   the GPU, or add IBM support to `solver_gpu` (offset the tridiagonal solve by `IBM_h`, build with
   `-D_USE_SOLVER_GPU`, and call it from `calculate_new_pressure`). Also check with `nsys`:
   arrays that are never prefetched (`mu_t`, `rhocp`, …) and per-step managed allocations
   (`outflow_ind` in `fillps`).
6. Making the IBM code decomposition-aware (offset `IBM_h` by `ijk_start(3)`) would allow splitting z.
7. **1-GPU runs blow up** (dims_in = 1 1): pressure checksum 3.8e14 in step 2, NaN by step 4, with both
   solvers (identical numbers), while 4 GPUs run fine. Something breaks when one rank is both
   inflow (SEM) and outflow, or exchanges with itself (`Receiver: 0 senders: 0`). Not investigated yet.
