# number of parallel compile jobs
NJOBS=8
# pressure solver on the GPU (1) or on the CPU (0, validated default)
USE_SOLVER_GPU=${USE_SOLVER_GPU:-0}

module load profile/candidate
module load nvhpc/25.3
module load hpcx-mpi/2.19
make clean APP=two_phase_cdi_marangoni && make -j $NJOBS USE_SOLVER_GPU=$USE_SOLVER_GPU ARCH=generic-gpu APP=two_phase_cdi_marangoni DO_DBG=0>& log_comp.o
tail -8 log_comp.o