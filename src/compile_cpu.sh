# CPU build (GNU) on Leonardo
# start from a clean environment: if nvhpc/hpcx-mpi are still loaded, mpif90 is the NVIDIA wrapper
module purge
module load profile/candidate
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2
# FFTW built with gcc-12.2.0 (same library used by the GPU build); overrides FFTW_HOME in targets/target.generic-gnu
FFTW_HOME=/leonardo/prod/spack/06/install/0.22/linux-rhel8-icelake/gcc-12.2.0/fftw-3.3.10-f3ms2lmqzjugb43xzrxjkvmvohumyruy
make clean APP=two_phase_cdi_marangoni && make ARCH=generic-gnu APP=two_phase_cdi_marangoni DO_DBG=0 FFTW_HOME=$FFTW_HOME >& log_comp.o
tail -8 log_comp.o
