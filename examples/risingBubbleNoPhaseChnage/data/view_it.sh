sleep 2
source ~/FluTAS/run/env/bin/activate
python3 write_xdmf.py
deactivate
paraview viewfld_DNS.xmf
