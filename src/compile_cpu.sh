module load gcc
module load openmpi
module load fftw
make clean APP=two_phase_cdi_marangoni && make ARCH=generic-gnu APP=two_phase_cdi_marangoni DO_DBG=0>& log_comp.o
tail -8 log_comp.o
