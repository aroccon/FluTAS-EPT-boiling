import csv
import sys
import glob
import numpy as np
import copy
import math
import os
import re
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sympy.solvers import solve,nsolve
from scipy.special import erf,erfc
from sympy import Symbol,Eq
import matplotlib.cm as cm
import sympy
import scipy.integrate as integrate
from scipy.optimize import minimize
params = {'legend.fontsize': 15,
          'legend.loc':'best',
          'figure.figsize': (15, 5),
          'lines.markerfacecolor':'none',
          'axes.labelsize': 15,
          'axes.titlesize': 15,
          'xtick.labelsize':15,
          'ytick.labelsize':15,
          'grid.alpha':0.6}
pylab.rcParams.update(params)
import vtk
from vtk.util import numpy_support as VN
from vtk.util.numpy_support  import vtk_to_numpy,numpy_to_vtk
from vtk.numpy_interface import algorithms as algs
from vtk.numpy_interface import dataset_adapter as dsa
from vtk import vtkStructuredPointsReader,vtkDataSetReader
from vtkmodules.vtkIOXML import (
    vtkXMLImageDataReader,
    vtkXMLImageDataWriter)
    
### Some subroutines

def read_and_write_film_data(folder,start):
    # first check if the .dat file exist
    if os.path.exists('./' + folder + '/data.dat'):
        print('The data file:  '+ './' + folder + '/data.dat' + ' does exist already: reading data file ...')
        times            = list()
        Nusselts          = list()
        input = open('./' + folder  + '/data.dat', 'r')
        lines = input.readlines()
        for line in lines[1:]:
            numbers = line.split()
            times.append(float(numbers[0]))
            Nusselts.append(float(numbers[1]))
        times = np.array(times)
        Nusselts = np.array(Nusselts)
    return times,Nusselts


def read_vti_scalar_data(file,field):
    reader = vtkXMLImageDataReader()
    reader.SetFileName(file)
    reader.Update()
    data = reader.GetOutput()
    ENT = VN.vtk_to_numpy(data.GetPointData().GetArray(field))
    dummy, pxdir,pydir = data.GetDimensions()
    x_0,y_0,z_0 = data.GetOrigin()
    spacing = data.GetSpacing()[0]
    ENT_grid = list()
    for i in range(0,pydir):
        ENT_grid.append(ENT[(i*pxdir):(i*pxdir+pxdir)])
    ENT_grid = np.array(ENT_grid)
    yp = np.linspace(y_0,y_0+(pxdir-1)*spacing,pxdir)
    zp = np.linspace(z_0,z_0+(pydir-1)*spacing,pydir)
    y_grid,z_grid = np.meshgrid(yp,zp)
    return y_grid,z_grid,ENT_grid
    
    

#The simulation parameters (https://doi.org/10.1080/10407790.2014.915681):
lamda_1 = 1
lamda_2 = lamda_1
cp1 = 200
cp2 = cp1
T_sat = 10
latent_heat = 1.0*1e4
T_wall = T_sat + 1
T_far = 10
rho1 = 5
rho2 = 200
g_acc = 9.81
mu1 = 0.005
mu2 = 0.1
alpha1 = lamda_1/(rho1*cp1)
alpha2 = lamda_2/(rho2*cp2)
sigma = 0.1

#some derived numbers and lengthscales

St = cp1*(T_wall-T_far)/latent_heat #stefan number
l_ref = 2*np.pi*np.sqrt(3*sigma/((rho2-rho1)*g_acc)) #unstable Taylor wave length
l_ast = np.sqrt(sigma/((rho2-rho1)*g_acc))  #The characteristic length scale

l_ast,l_ref,St


# Plot some VTK files

folder        = 'data'
files_to_plot = ['000050.vti','000075.vti','000100.vti']
color         = 'hotpink'

fig, ax = plt.subplots(1,3,figsize=(9, 4),sharey=True)
raster = False
x_loc = 0.001
y_loc = 0.073
offset = l_ref
axis = 1
ax[0].set_ylabel('y [m]')

for a in ax:
    a.set_aspect(1.0)
    a.ticklabel_format(axis = 'both',style = 'sci',useMathText=True, scilimits=(0,0))
    a.set_xlim([0,offset])
    a.set_ylim([0,0.08])
    a.set_xlabel('x [m]')

for file in reversed(files_to_plot):
    x,y,phi = read_vti_scalar_data(folder + '/' + file ,'psi')
    x_ = np.concatenate((x,np.flip(-x+offset, axis=axis)),axis=axis)
    y_ = np.concatenate((y,np.flip(y, axis=axis)),axis=axis)
    phi_ = np.concatenate((phi,np.flip(phi, axis=axis)),axis=axis)
    ax[files_to_plot.index(file)].contour(x_,y_,phi_,levels=[0.5],linestyles=['solid'],colors=[color])
    time = float(file.split('.')[0].lstrip('0'))/100 # in ms
    time_text = r't= {:0.2f} s'.format(time)
    ax[files_to_plot.index(file)].text(x_loc,y_loc,time_text,fontsize=12)

ax[2].plot([1,1],[1,1],c=color,linestyle='solid', label = r'$\phi$ = 0.5')

lgd = ax[2].legend(loc='upper center', bbox_to_anchor=(1.4, 0.7),
          fancybox=False, shadow=False, ncol=1)

fig.tight_layout()
fig.savefig('BubbleDetachment.pdf', format='pdf', bbox_inches='tight')


# Plot the Nusselt numbers

times,Nusselts = read_and_write_film_data(folder,0.0)


# Nussel number correlation for comparison is taken from:  https://doi.org/10.1016/0017-9310(81)90094-6
# See also here https://doi.org/10.1080/10407790190053824 and https://doi.org/10.1080/10407782.2011.561079

#Guo 2011
Gr = (rho2/rho1-1)*g_acc*l_ast**3*rho1**2/mu1**2  #Grashof number  
beta = 1/(latent_heat/(cp1*(T_wall-T_sat)))  #For Klimenko corr
f1 = 0.89*(beta)**(-1/3)   #For Klimenko corr -> because smaller than 0.71
Pr = (mu1/rho1)/alpha1  # Prantl number
Nu_Kli = 0.19*(Gr)**(1/3)*(Pr**(1/3))*f1 # Klimenko corr


fig, ax = plt.subplots(figsize=(5, 4))
resolution = 500
skip = int(len(times)/resolution)


ax.plot(times[::skip],Nusselts[::skip],c = '0',linestyle='solid', label = r'Simulation')
ax.plot([times.min(),times.max()],[Nu_Kli,Nu_Kli],c='hotpink',label='Klimenko (1981)')


ax.legend(fancybox=False, shadow=False)

ax.set_ylim([0,4])
ax.set_xlim([times.min(),times.max()])
ax.ticklabel_format(axis = 'both',style = 'sci',useMathText=True, scilimits=(0,0))
ax.set_xlabel('Time [s]')
ax.set_ylabel(r'Nusselt [-]')
fig.savefig('NusseltNumbers.pdf', format='pdf', bbox_inches='tight',transparent=True)


