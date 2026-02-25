module load profile/candidate
module load nvhpc/25.3
module load hpcx-mpi/2.19
make clean APP=two_phase_cdi_marangoni && make ARCH=generic-gpu APP=two_phase_cdi_marangoni DO_DBG=0>& log_comp.o
tail -8 log_comp.o