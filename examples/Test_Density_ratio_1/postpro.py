import os
import numpy as np
import sys
import glob
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {'legend.fontsize': 15,
          'legend.loc':'best',
          'figure.figsize': (15, 5),
         'lines.markerfacecolor':'none',
         'axes.labelsize': 17,
         'axes.titlesize': 17,
         'xtick.labelsize':15,
         'ytick.labelsize':15,
         'grid.alpha':0.6}
pylab.rcParams.update(params)
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sympy.solvers import solve,nsolve
from sympy import Symbol,Eq
import sympy
import re
from scipy.special import erf,erfc

def read_single_field_binary(filenamei,iskip,folder):
    #
    # setting up some parameters
    #
    iprecision = 8            # precision of the real-valued data
    r0 = np.array([0.,0.,0.]) # domain origin
    non_uniform_grid = True
    precision  = 'float64'
    if(iprecision == 4): precision = 'float32'
    #
    # read geometry file
    #
    geofile  = folder + "geometry.out"
    geo = np.loadtxt(geofile, comments = "!", max_rows = 2)
    ng = geo[0,:].astype('int')
    l  = geo[1,:]
    dl = l/(1.*ng)
    #
    # read and generate grid
    #
    xp = np.arange(r0[0]+dl[0]/2.,r0[0]+l[0],dl[0]) # centered  x grid
    yp = np.arange(r0[1]+dl[1]/2.,r0[1]+l[1],dl[1]) # centered  y grid
    zp = np.arange(r0[2]+dl[2]/2.,r0[2]+l[2],dl[2]) # centered  z grid
    xu = xp + dl[0]/2.                              # staggered x grid
    yv = yp + dl[1]/2.                              # staggered y grid
    zw = zp + dl[2]/2.                              # staggered z grid
    if(non_uniform_grid):
        f   = open(folder + 'grid.bin','rb')
        grid_z = np.fromfile(f,dtype=precision)
        f.close()
        grid_z = np.reshape(grid_z,(ng[2],4),order='F')
        zp = r0[2] + np.transpose(grid_z[:,2]) # centered  z grid
        zw = r0[2] + np.transpose(grid_z[:,3]) # staggered z grid
    #
    # read binary file
    #
    n           = (ng[:]/iskip[:]).astype(int)
    data        = np.zeros([n[0],n[1],n[2]])
    fld         = np.fromfile(filenamei,dtype=precision)
    data[:,:,:] = np.reshape(fld,(n[0],n[1],n[2]),order='F')
    #
    # reshape grid
    #
    xp = xp[0:ng[0]:iskip[0]]
    yp = yp[0:ng[1]:iskip[1]]
    zp = zp[0:ng[2]:iskip[2]]
    xu = xu[0:ng[0]:iskip[0]]
    yv = yv[0:ng[1]:iskip[1]]
    zw = zw[0:ng[2]:iskip[2]]
    return data,xp,yp,zp,xu,yv,zw

def get_avg_psi_loc(psi,xp,yp,zp,target):
    loc = 0.0
    count = 0
    psi = np.array(psi)    
    for i in range(0,len(xp)):
        for j in range(0,len(yp)):
            for k in range(0,len(zp)):
                if(psi[i,j,k]<=target):
                    break
            if((k==len(zp)-1)or(k==0)):
                loc = loc + 0.0
                count = count + 1
            else:
                alpha = (target-psi[i,j,k-1])/(psi[i,j,k]-psi[i,j,k-1])
                loc = loc + zp[k]*alpha + zp[k-1]*(1.0-alpha)
                count = count + 1
    loc = loc/count
    return loc

def get_analytical_sol():
    lamda_1 = 0.0005
    lamda_2 = 0.005
    cp1 = 200.0
    cp2 = 700.0
    T_sat = 0.0
    latent_heat = 10000.0
    T_wall = T_sat
    T_far = 10.0
    rho1 = 0.001
    rho2 = 0.001
    mu1 = 0.00001
    mu2 = 0.001
    alpha1 = lamda_1/(rho1*cp1)
    alpha2 = lamda_2/(rho2*cp2)
    A = cp1*(T_wall-T_sat)/(latent_heat*sympy.sqrt(sympy.pi))
    eps = Symbol('eps', real=True)
    B = ((T_far-T_sat)*cp1*lamda_2*sympy.sqrt(alpha1)*sympy.exp(-(eps**2)*((alpha1*rho1**2)/(alpha2*rho2**2))))
    C = (latent_heat*lamda_1*sympy.sqrt(sympy.pi*alpha2)*sympy.erfc(eps*((sympy.sqrt(alpha1)*rho1)/(sympy.sqrt(alpha2)*rho2))))
    equation1 = A - sympy.exp(eps**2)*sympy.erf(eps)*(eps-B/C)
    eps_num = float(nsolve(equation1, eps ,1, prec=25))
    time = np.linspace(0,4,1000)
    front = 2*eps_num*np.sqrt(lamda_1*time/(rho1*cp1))
    return time,front

def get_analytical_pos(time):
    lamda_1 = 0.0005
    lamda_2 = 0.005
    cp1 = 200.0
    cp2 = 700.0
    T_sat = 0.0
    latent_heat = 10000.0
    T_wall = T_sat
    T_far = 10.0
    rho1 = 0.001
    rho2 = 0.001
    mu1 = 0.00001
    mu2 = 0.001
    alpha1 = lamda_1/(rho1*cp1)
    alpha2 = lamda_2/(rho2*cp2)
    A = cp1*(T_wall-T_sat)/(latent_heat*sympy.sqrt(sympy.pi))
    eps = Symbol('eps', real=True)
    B = ((T_far-T_sat)*cp1*lamda_2*sympy.sqrt(alpha1)*sympy.exp(-(eps**2)*((alpha1*rho1**2)/(alpha2*rho2**2))))
    C = (latent_heat*lamda_1*sympy.sqrt(sympy.pi*alpha2)*sympy.erfc(eps*((sympy.sqrt(alpha1)*rho1)/(sympy.sqrt(alpha2)*rho2))))
    equation1 = A - sympy.exp(eps**2)*sympy.erf(eps)*(eps-B/C)
    eps_num = float(nsolve(equation1, eps ,1, prec=25))
    front = 2*eps_num*np.sqrt(lamda_1*time/(rho1*cp1))
    return front

def get_analytical_temps(time):
    lamda_1 = 0.0005
    lamda_2 = 0.005
    cp1 = 200.0
    cp2 = 700.0
    T_sat = 0.0
    latent_heat = 10000.0
    T_wall = T_sat
    T_far = 10.0
    rho1 = 0.001
    rho2 = 0.001
    mu1 = 0.00001
    mu2 = 0.001
    alpha1 = lamda_1/(rho1*cp1)
    alpha2 = lamda_2/(rho2*cp2)
    A = cp1*(T_wall-T_sat)/(latent_heat*sympy.sqrt(sympy.pi))
    eps = Symbol('eps', real=True)
    B = ((T_far-T_sat)*cp1*lamda_2*sympy.sqrt(alpha1)*sympy.exp(-(eps**2)*((alpha1*rho1**2)/(alpha2*rho2**2))))
    C = (latent_heat*lamda_1*sympy.sqrt(sympy.pi*alpha2)*sympy.erfc(eps*((sympy.sqrt(alpha1)*rho1)/(sympy.sqrt(alpha2)*rho2))))
    equation1 = A - sympy.exp(eps**2)*sympy.erf(eps)*(eps-B/C)
    eps_num = float(nsolve(equation1, eps ,1, prec=25))
    loc = 2*eps_num*np.sqrt(lamda_1*time/(rho1*cp1))
    x1 = np.linspace(0,loc,100)
    T1 = T_wall+((T_sat-T_wall)/erf(eps_num))*erf(x1/(2*np.sqrt(alpha1*time)))
    x2 = np.linspace(loc,0.8,400)
    T2 = T_far-((T_far-T_wall)/erfc(eps_num*((np.sqrt(alpha1)*rho1)/(np.sqrt(alpha2)*rho2))))*erfc(x2/(2*np.sqrt(alpha2*time))+eps_num*(rho1-rho2)/rho2*np.sqrt(alpha1/alpha2))
    return x1,T1,x2,T2



def get_time(file):
    stamp = ' ' + re.findall(r'\d+',file)[-1].lstrip("0") + ' '
    time = 0
    found = False
    with open(r"log.o", 'r') as fp:
        for l_no, line in enumerate(fp):
            # search string
            if ((stamp in line)and(r'Timestep #' in line)):
                found = True
                time = float(line.split()[-1])
                break
    if not found:
        print('No time was found for timestep: '+stamp)
    return time

fig, axs = plt.subplots(1,4,figsize=(20, 5))

i = 0
start = 0.05
sub_dir = 'data'
print('simulation starts at time [s]:')
print(start)

fig, ax = plt.subplots(figsize=(8, 5))
divs_ibm = list()
divs_u = list()
divs_l = list()
times_ibm = list()

folder = './' + sub_dir + '/'
files = glob.glob('./' + sub_dir + '/psi_fld_*.bin')
files.sort()
i = 0
for file in files:
    time = get_time(file)+start
    filenamei   = file
    iskip       = np.array([1,1,1]).astype(int)
    div,xp,yp,zp,xu,yv,zw = read_single_field_binary(filenamei,iskip,folder)
    times_ibm.append(time)
    divs_ibm.append(get_avg_psi_loc(div,xp,yp,zp,0.5))
    divs_l.append(get_avg_psi_loc(div,xp,yp,zp,0.1))
    divs_u.append(get_avg_psi_loc(div,xp,yp,zp,0.9))
    i = i + 1

ax.plot(times_ibm,divs_ibm,c = 'r',label = r'CDI ($\phi = 0.5$)',linestyle = 'dashed')
print('time points:')
print(times_ibm)
print('location points:')
print(divs_ibm)
ax.fill_between(times_ibm,divs_u,divs_l,color = '0.9',label = 'CDI ($0.1 < \phi < 0.9$)')
time,front = get_analytical_sol()
ax.plot(time,front,c = '0',label = 'Analytical',linestyle = 'solid')
ax.legend()
ax.set_xlabel(r'time [s]')
ax.set_ylabel(r'interface position [m]')
#ax.set_xlim([0,0.08e-3])
#ax.set_ylim([0.0,0.0001])
ax.ticklabel_format(axis = 'both',style = 'sci',useMathText=True, scilimits=(0,0))
fig.tight_layout()
fig.savefig('loc_vs_t.pdf', format='pdf')#,dpi = 600)

field_nums = [0,25,50,75,100]#,50]
fig, ax = plt.subplots(figsize=(8, 5))

for file_num in field_nums:
    folder = './' + sub_dir + '/'
    files = glob.glob('./' + sub_dir + '/tmp_fld_*.bin')
    files.sort()
    filenamei   = files[file_num]
    time = get_time(filenamei)+start
    pos = get_analytical_pos(time)
    x1,T1,x2,T2 = get_analytical_temps(time)
    if(file_num==0):
        ax.plot(x1,T1,c = '0',label = 'Analytical',linestyle = 'solid')#,zorder = 4)
        ax.plot(x2,T2,c = '0',linestyle = 'solid')#,zorder = 4)
    else:
        ax.plot(x1,T1,c = '0',linestyle = 'solid')#,zorder = 4)
        ax.plot(x2,T2,c = '0',linestyle = 'solid')#,zorder = 4)
    #ax.set_xlim([0,0.0001])    
    iskip       = np.array([1,1,1]).astype(int)
    div,xp,yp,zp,xu,yv,zw = read_single_field_binary(filenamei,iskip,folder)
    div = np.array(div)
    div = div.mean(axis=(0, 1))
    div = div
    ax.plot(zp,div,c = 'r',label = 'time {:0.2f}'.format(time),linestyle = 'dashed')#,zorder = 5)   
    #ax.set_xticks([0,0.005,0.01,0.015,0.02])
    #ax.set_yticks([0,0.005,0.01,0.015,0.02])
    ax.ticklabel_format(axis = 'both',style = 'sci',useMathText=True, scilimits=(0,0))

ax.set_xlabel(r'Walldistance [m]')
ax.set_ylabel(r'$T-T_{sat}$ [K]')
ax.legend()
fig.tight_layout()
fig.savefig('temps.pdf', format='pdf') #,dpi = 600)
