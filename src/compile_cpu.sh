# module load gcc
#module load openmpi
#module load fftw
make clean APP=two_phase_acdi_boiling && make ARCH=generic-gnu APP=two_phase_acdi_boiling DO_DBG=0>& log_comp.o
tail -8 log_comp.o
