module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2 
module load fftw/3.3.10--openmpi--4.1.6--gcc--12.2.0-spack0.22
nmake clean APP=two_phase_cdi_marangoni && make ARCH=generic-gnu APP=two_phase_cdi_marangoni DO_DBG=0>& log_comp.o
tail -8 log_comp.o
