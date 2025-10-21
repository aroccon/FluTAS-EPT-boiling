import os
import numpy as np
import sys
import glob

rho1 = float(input('Rho1:'))
rho2 = float(input('Rho2:'))
cp1 = float(input('Cp1:'))
cp2 = float(input('Cp2:'))
tmax = float(input('T_max:'))
hlat = float(input('lat_heat:'))

h02 = -(tmax*(rho1*cp1-rho2*cp2)+rho1*hlat)/(rho1-rho2)
h01 = h02 + hlat
print('   h01      h02   ')
print('{:0.2f}'.format(h01)+'     '+'{:0.2f}'.format(h02))
