import csv
import sys
import glob
import numpy as np
import copy
import os
import re
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sympy.solvers import solve,nsolve
from scipy.special import erf,erfc
from sympy import Symbol,Eq
import sympy
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

sys.path.insert(0,'./2014_Adelsberger_reference_data')

import plot_reference_data

################################
#   required subroutines       #
################################


def read_single_field_binary(filenamei,iskip,folder):
    # Reading in binary data
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
    

def get_time(file,folder):
    # read time from logfile or outputfile recoder (located in ./data/ ) 
    filename = file.split('/')[-1].split('\\')[-1]
    stamp = ' ' + re.findall(r'\d+',file)[-1].lstrip("0") + ' '
    time = 0
    found = False
    with open(folder+"log_visu_3d.out", 'r') as fp:
        for l_no, line in enumerate(fp):
            # search string
            if (filename in line):
                found = True
                time = float(line.split()[-2])
                break
    if not found:
        print('No time was found for timestep: '+stamp+' in file log_visu_3d.out...trying log.out')
        with open(folder+"log.out", 'r') as fp:
            for l_no, line in enumerate(fp):
            # search string
                if ((stamp in line)and(r'Timestep #' in line)):
                    found = True
                    time = float(line.split()[-1])
                    break
    if not found:
        print('No time was found for timestep: '+stamp)
    return time

def get_center_of_mass(psi,zp,cutoff):
    # get center of mass of the whole bubble
    sum_1 = 0
    sum_2 = 0
    for i in range(0,len(psi)):
        for j in range(0,len(psi[i])):
            for k in range(0,len(psi[i][j])):
                if(psi[i][j][k]>cutoff):
                    sum_1 = sum_1 + psi[i][j][k]*zp[k]
                    sum_2 = sum_2 + psi[i][j][k]
    center = sum_1/sum_2
    return center

def get_riseVelocity(psi,vez,cutoff):
    # get integral bubble velocity data
    sum_1 = 0
    sum_2 = 0
    for i in range(0,len(vez)):
        for j in range(0,len(vez[i])):
            for k in range(0,len(vez[i][j])-1):
                if((0.5*(psi[i][j][k+1]+psi[i][j][k]))>cutoff):
                    sum_1 = sum_1 + (0.5*(psi[i][j][k+1]+psi[i][j][k]))*vez[i][j][k]
                    sum_2 = sum_2 + (0.5*(psi[i][j][k+1]+psi[i][j][k]))
    velocity = sum_1/sum_2
    return velocity

def read_sphericity(sub_dir,t_ref):
    #read bubble sphericity from tagging file
    dat_sphereicity = list()
    t_sphereicity = list()
    bubb_data = list()
    timesteps = set()
    input = open('./'+sub_dir+'/post/concatenate/concatenated_bubbles.dat', 'r')
    lines = input.readlines()
    data_names = lines[0].split()
    for line in lines[1:]:
        bubble = list()
        split_line = line.split()
        timesteps.add(int(split_line[0]))
        for entry in split_line:
            bubble.append(float(entry))
        bubb_data.append(bubble)
    input.close()
    for bubble in bubb_data:
      if(bubble[data_names.index('bubble_ID')]==1.0):
          t_sphereicity.append(bubble[data_names.index('Time')])
          dat_sphereicity.append(bubble[data_names.index('surf')])
    t_sphereicity = np.array(t_sphereicity)
    dat_sphereicity = np.array(dat_sphereicity)
    t_sphereicity = t_sphereicity/t_ref
    dat_sphereicity = dat_sphereicity[0]/dat_sphereicity
    return t_sphereicity,dat_sphereicity

def read_bubble_center_and_velocity(sub_dir,t_ref,u_ref):
    #read bubble center and bubble velocity from binary file
    dat_centerOfMass = list()
    t_centerOfMass = list()
    dat_riseVelocity = list()
    t_riseVelocity = list()
    folder = './' + sub_dir + '/'
    files = glob.glob('./' + sub_dir + '/psi_fld_*.bin')
    files.sort()
    files_vel = glob.glob('./' + sub_dir + '/vez_fld_*.bin')
    files_vel.sort()
    i = 0
    for i in range(0,len(files_vel)):
        iskip       = np.array([1,1,1]).astype(int)
        time = get_time(files_vel[i],folder)
        t_riseVelocity.append(time)
        t_centerOfMass.append(time)
        print('Reading '+files[i]+' (',(i+1),'/',len(files),') : ', time, '...')
        filenamei   = files[i]
        psi,xp,yp,zp,xu,yv,zw = read_single_field_binary(filenamei,iskip,folder)
        centerOfMass = get_center_of_mass(psi,zp,0.01)
        dat_centerOfMass.append(centerOfMass)
        print('done.')
        filenamei_v = files_vel[i]
        print('Reading '+filenamei_v+' (',(i+1),'/',len(files_vel),') : ', time, '...')
        vez,xp,yp,zp,xu,yv,zw = read_single_field_binary(filenamei_v,iskip,folder)
        riseVelocity = get_riseVelocity(psi,vez,0.01)
        dat_riseVelocity.append(riseVelocity)
        print('done.')
        i = i + 1
    t_centerOfMass = np.array(t_centerOfMass)
    dat_centerOfMass = np.array(dat_centerOfMass)
    t_centerOfMass = t_centerOfMass/t_ref
    t_riseVelocity = np.array(t_riseVelocity)
    dat_riseVelocity = np.array(dat_riseVelocity)
    t_riseVelocity = t_riseVelocity/t_ref
    dat_riseVelocity = dat_riseVelocity/u_ref
    return t_centerOfMass,dat_centerOfMass,t_riseVelocity,dat_riseVelocity


def read_bubble_center_and_velocity2(sub_dir,t_ref,u_ref): 
    #read bubble center and bubble velocity from tagging file
    dat_centerOfMass = list()
    t_centerOfMass = list()
    dat_riseVelocity = list()
    t_riseVelocity = list()
    bubb_data = list()
    timesteps = set()
    input = open('./'+ sub_dir +'/post/concatenate/concatenated_bubbles.dat', 'r')
    lines = input.readlines()
    data_names = lines[0].split()
    for line in lines[1:]:
        bubble = list()
        split_line = line.split()
        timesteps.add(int(split_line[0]))
        for entry in split_line:
            bubble.append(float(entry))
        bubb_data.append(bubble)
    input.close()
    for bubble in bubb_data:
      if(bubble[data_names.index('bubble_ID')]==1.0):
          t_centerOfMass.append(bubble[data_names.index('Time')])
          t_riseVelocity.append(bubble[data_names.index('Time')])
          dat_centerOfMass.append(bubble[data_names.index('zpos')])
          dat_riseVelocity.append(bubble[data_names.index('wvel')])
    t_centerOfMass = np.array(t_centerOfMass)
    dat_centerOfMass = np.array(dat_centerOfMass)
    t_centerOfMass = t_centerOfMass/t_ref
    t_riseVelocity = np.array(t_riseVelocity)
    dat_riseVelocity = np.array(dat_riseVelocity)
    t_riseVelocity = t_riseVelocity/t_ref
    dat_riseVelocity = dat_riseVelocity/u_ref
    return t_centerOfMass,dat_centerOfMass,t_riseVelocity,dat_riseVelocity


def read_and_write_data(sub_dir,t_ref,u_ref):
    # first check if the .dat file exists
    if os.path.exists('./' + sub_dir + '/'+sub_dir+'.dat'):
        print('The data file:  '+ './' + sub_dir + '/data.dat' + ' does exist already: reading data file ...')
        t_centerOfMass = list()
        dat_centerOfMass = list()
        t_riseVelocity = list()
        dat_riseVelocity = list()
        t_sphereicity = list()
        dat_sphereicity = list()
        input = open('./' + sub_dir + '/'+sub_dir+'.dat', 'r')
        lines = input.readlines()
        for line in lines[1:]:
            numbers = line.split()
            t_centerOfMass.append(float(numbers[0]))
            dat_centerOfMass.append(float(numbers[1]))
            t_riseVelocity.append(float(numbers[0]))
            dat_riseVelocity.append(float(numbers[2]))
            t_sphereicity.append(float(numbers[0]))
            dat_sphereicity.append(float(numbers[3]))
        t_centerOfMass = np.array(t_centerOfMass)
        dat_centerOfMass = np.array(dat_centerOfMass)
        t_riseVelocity = np.array(t_riseVelocity)
        dat_riseVelocity = np.array(dat_riseVelocity)
        t_sphereicity = np.array(t_sphereicity)
        dat_sphereicity = np.array(dat_sphereicity)
    else:  # if it doesn't exist read data files
        print('The data file:  '+ './' + sub_dir + '/data.dat' + ' does not exist yet: reading binary files ...')
        t_centerOfMass,dat_centerOfMass,t_riseVelocity,dat_riseVelocity = read_bubble_center_and_velocity2(sub_dir,t_ref,u_ref)
        t_sphereicity,dat_sphereicity = read_sphericity(sub_dir,t_ref)
        with open('./' + sub_dir + '/data.dat', 'w') as file:
            file.write('Time     CenterOfMass      BubbleVelocity     Sphericity\n')
            for i in range(0,len(t_centerOfMass)):
                file.write(str(t_sphereicity[i])+'     '+str(dat_centerOfMass[i])+'     '+str(dat_riseVelocity[i])+'     '+str(dat_sphereicity[i]))
                if i!=(len(t_centerOfMass)-1): file.write('\n')
    return t_centerOfMass,dat_centerOfMass,t_riseVelocity,dat_riseVelocity,t_sphereicity,dat_sphereicity

################################
#     Simulation setting       #
################################

#non-dimensional data
Re = 35.0
Fr = 1.0
We = 1.0
rho_rat = 10.0
mu_rat = 10.0

#dimensions from the simulation

#dimensional
rho1 = 100.0
g = -0.98
d0 = 0.5


#derived
rho2 = rho1*rho_rat
lx = 2.0*d0
ly = lx
lz = 4.0*d0
l_ref = d0
u_ref = np.sqrt(abs(g)*d0)
t_ref = np.sqrt(d0/abs(g))
mu1 = rho1*u_ref*l_ref/Re
mu2 = mu1*mu_rat
sigma = (rho1*l_ref*u_ref**2)/We

################################
#         Visualize            #
################################


fig, axs = plt.subplots(1,3,figsize=(22, 5))

axs = plot_reference_data.plotAdelsberger(axs)

entries = [['data','hotpink','dashed','Present']]

for entry in entries:
    BC_t, BC , RV_t , RV , SPH_t , SPH = read_and_write_data(entry[0],t_ref,u_ref)
    axs[0].plot(BC_t,BC,c=entry[1],linestyle=entry[2])
    axs[1].plot(SPH_t,SPH,c=entry[1],linestyle=entry[2])
    axs[2].plot(RV_t,RV,c=entry[1],linestyle=entry[2],label=entry[3])


axs[0].set_xlabel('Time [-]')
axs[1].set_xlabel('Time [-]')
axs[2].set_xlabel('Time [-]')
axs[0].set_ylabel('Center of mass [-]')
axs[1].set_ylabel('Sphericity [-]')
axs[2].set_ylabel('Bubble velocity [-]')
axs[2].legend()


fig.suptitle('Bubble rise simulation compared to  Adelsberger (2014)',fontsize=20)
fig.savefig('Bubble_rise_comparison.pdf', format='pdf', bbox_inches='tight')




