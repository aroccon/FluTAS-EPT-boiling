. ~/add_nvidia.bash
make clean APP=two_phase_acdi_boiling && make ARCH=generic-gpu APP=two_phase_acdi_boiling DO_DBG=0>& log_comp.o
tail -8 log_comp.o
