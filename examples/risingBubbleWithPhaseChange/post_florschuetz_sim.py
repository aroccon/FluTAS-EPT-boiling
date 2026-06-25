import csv
import sys
import glob
import numpy as np
import copy
import os
import re
import math
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sympy.solvers import solve,nsolve
from scipy.special import erf,erfc
from sympy import Symbol,Eq
import sympy
import scipy.integrate as integrate
from scipy.optimize import minimize
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
import vtk
from vtk.util import numpy_support as VN
from vtk.util.numpy_support  import vtk_to_numpy,numpy_to_vtk
from vtk.numpy_interface import algorithms as algs
from vtk.numpy_interface import dataset_adapter as dsa
from vtk import vtkStructuredPointsReader,vtkDataSetReader
from vtkmodules.vtkIOXML import (
    vtkXMLImageDataReader,
    vtkXMLImageDataWriter)


#some subroutines

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

def integrand(zeta, eps, rho1, rho2):
    return np.exp(-(eps**2.0)*((1.0-zeta)**(-2.0)-2.0*(1.0-rho1/rho2)*zeta-1.0))

def opti(eps,rho1,rho2,cp1,cp2,T_far,T_sat,latent_heat):
    A = cp2*rho2*(T_far-T_sat)/(rho1*(latent_heat+(cp2-cp1)*(T_far-T_sat)))
    C = integrate.quad(integrand,0,1,args = (eps,rho1,rho2),epsabs = 1e-20)[0]
    D = 2.0*C*eps**2
    error = abs(A - D)
    return error

def get_time(file,folder):
    filename = file.split('/')[-1].split('\\')[-1]
    stamp = ' ' + re.findall(r'\d+',file)[-1].lstrip("0") + ' '
    time = 0
    found = False
    if(stamp==''):
        time = 0.0
        found = True
    if not found:
        with open(folder+"log.out", 'r') as fp:
            for l_no, line in enumerate(fp):
                # search string 
                if ((stamp in line)and(r'Timestep #' in line)):
                    found = True
                    time = float(line.split()[-1])
                    break
    if not found:
        print('Warning: No time was found for timestep: '+stamp+' in file log.out...trying log_visu_3d.out')
        if os.path.exists(folder+"log_visu_3d.out"):
            with open(folder+"log_visu_3d.out", 'r') as fp:
                for l_no, line in enumerate(fp):
                    # search string
                    if (filename in line):
                        found = True
                        time = float(line.split()[-2])
                        break
    if not found:
        print('No time was found for timestep: '+stamp)
    return time

def get_inflow_data(folder):
    output = [[0.0],[0.0],[0.0],[0.0]]
    if not os.path.exists(folder+"log.out"):
        sys.exit('ERROR: No file '+folder+"log.out"+' was not found...exiting')
    with open(folder+"log.out", 'r') as fp:
        for l_no, line in enumerate(fp):
            # search string
            if (r'Inflow velocity:' in line):
                dist_old = output[3][-1]
                dt       = float(line.split()[-1])
                time     = float(line.split()[-2])
                invel    = float(line.split()[-3])
                dist_new = dist_old - dt*invel
                output[0].append(time)
                output[1].append(dt)
                output[2].append(invel)
                output[3].append(dist_new)
    return output

def get_corr_vel(time,output):
    corr_vel = output[2][output[0].index(time)]
    return corr_vel

def get_corr_loc(time,output):
    corr_loc = output[3][output[0].index(time)]
    return corr_loc

def create_1D_lineout(psi,xp,yp,zp):
    loc     = 0.0
    maxdist = np.sqrt(xp[-1]**2 + yp[-1]**2 + zp[-1]**2)
    psi     = np.array(psi)
    nx      = int(np.sqrt(2)*len(xp))-1
    xs      = np.linspace(0.0,maxdist,nx)
    psi_int = np.zeros(nx)
    wgt_int = np.zeros(nx)
    dx      = (xs[-1]-xs[0])/(len(xs)-1)
    for i in range(0,len(xp)):
        for j in range(0,len(yp)):
            for k in range(0,len(zp)):
                dist               = np.sqrt(xp[i]**2 + yp[j]**2 + zp[k]**2)
                if (dist == xs[-1]):
                    im_int         = len(xs) - 2
                else:
                    im_int         = math.floor(dist/dx)
                alpha              = (dist-xs[im_int])/(xs[im_int+1]-xs[im_int])
                psi_int[im_int]   += (1.0-alpha)*psi[i,j,k]
                psi_int[im_int+1] += (    alpha)*psi[i,j,k]
                wgt_int[im_int]   += (1.0-alpha)
                wgt_int[im_int+1] += (    alpha)
    for i in range(0,len(psi_int)):
        if((psi_int[i] == 0.0)&(wgt_int[i] == 0.0)):
            wgt_int[i] = 1.0
    psi_int = psi_int/wgt_int
    return xs,psi_int

def get_avg_psi_loc(xs,psi_int,target):
    loc     = 0.0
    for i in range(0,len(xs)):
        if(psi_int[i]<=target):
            break
    if((i==len(xs)-1)or(i==0)):
        loc = 0.0
    else:
        alpha = (target-psi_int[i-1])/(psi_int[i]-psi_int[i-1])
        loc = xs[i]*alpha + xs[i-1]*(1.0-alpha)
    return loc

def get_d_z(psi,zp):
    psi_center = psi[0][0]
    lower,upper = get_upper_and_lower_psi_loc(zp,psi_center,0.5)
    d_z = upper-lower
    return d_z

def get_d_r(psi,xp,yp,zp):
    d_r     = 0.0
    r_r_max = 0.0
    for k in range(0,len(xp)):
        z_slice    = psi[:,:,k]
        xs,psi_int = create_1D_lineout_from_slice(z_slice,xp,yp,zp)
        r_r        = get_avg_psi_loc(xs,psi_int,0.5)
        r_r_max = max(r_r_max,r_r)
    d_r = 2.0*r_r_max
    return d_r

def get_upper_and_lower_psi_loc(xs,psi_int,target):
    #first lower
    lower     = 0.0
    for i in range(0,len(xs)):
        if(psi_int[i]>=target):
            break
    if((i==len(xs)-1)or(i==0)):
        lower = 0.0
    else:
        alpha = (target-psi_int[i-1])/(psi_int[i]-psi_int[i-1])
        lower = xs[i]*alpha + xs[i-1]*(1.0-alpha)
    # now upper
    upper     = max(xs)
    for i in reversed(range(0,len(xs))):
        if(psi_int[i]>=target):
            break
    if((i==len(xs)-1)or(i==0)):
        upper = max(xs)
    else:
        alpha = (target-psi_int[i+1])/(psi_int[i]-psi_int[i+1])
        upper = xs[i]*alpha + xs[i+1]*(1.0-alpha)
    return lower,upper

def create_1D_lineout_from_slice(z_slice,xp,yp,zp):
    loc     = 0.0
    maxdist = np.sqrt(xp[-1]**2 + yp[-1]**2)
    z_slice = np.array(z_slice)
    nx      = int(np.sqrt(2)*len(xp))-1
    xs      = np.linspace(0.0,maxdist,nx)
    psi_int = np.zeros(nx)
    wgt_int = np.zeros(nx)
    dx      = (xs[-1]-xs[0])/(len(xs)-1)
    for i in range(0,len(xp)):
        for j in range(0,len(yp)):
            dist               = np.sqrt(xp[i]**2 + yp[j]**2)
            if (dist == xs[-1]):
                im_int         = len(xs) - 2
            else:
                im_int         = math.floor(dist/dx)
            alpha              = (dist-xs[im_int])/(xs[im_int+1]-xs[im_int])
            psi_int[im_int]   += (1.0-alpha)*z_slice[i,j]
            psi_int[im_int+1] += (    alpha)*z_slice[i,j]
            wgt_int[im_int]   += (1.0-alpha)
            wgt_int[im_int+1] += (    alpha)
    for i in range(0,len(psi_int)):
        if((psi_int[i] == 0.0)&(wgt_int[i] == 0.0)):
            wgt_int[i] = 1.0
    psi_int = psi_int/wgt_int
    return xs,psi_int

def get_center_of_mass(psi,zp,cutoff):
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

def read_bubble_data(sub_dir,corr_list,start):
    iskip       = np.array([1,1,1]).astype(int)
    t_               = list()
    dat_centerOfMass = list()
    dat_riseVelocity = list()
    dat_d_z          = list()
    dat_d_r          = list()
    folder = './' + sub_dir + '/'
    files = glob.glob('./' + sub_dir + '/psi_fld_*.bin')
    files.sort()
    files_vel = glob.glob('./' + sub_dir + '/vez_fld_*.bin')
    files_vel.sort()
    i = 0
    for i in range(0,len(files_vel)):
        time = get_time(files_vel[i],folder) #simulation time
        z_corr   = get_corr_loc(time,corr_list)
        vel_corr = get_corr_vel(time,corr_list)
        t_.append(time+start) #physical time
        print('Reading '+files[i]+' (',(i+1),'/',len(files),') : ', time, '...')
        filenamei   = files[i]
        psi,xp,yp,zp,xu,yv,zw = read_single_field_binary(filenamei,iskip,folder)
        zp = zp + z_corr
        print('Calculate center of mass...')
        centerOfMass = get_center_of_mass(psi,zp,0.01)
        dat_centerOfMass.append(centerOfMass)
        print('Calculate bubble size...')
        d_z = get_d_z(psi,zp)
        d_r = get_d_r(psi,xp,yp,zp)
        dat_d_z.append(d_z)
        dat_d_r.append(d_r)
        print('done.')
        filenamei_v = files_vel[i]
        print('Reading '+filenamei_v+' (',(i+1),'/',len(files_vel),') : ', time, '...')
        vez,xp,yp,zp,xu,yv,zw = read_single_field_binary(filenamei_v,iskip,folder)
        vez = vez - vel_corr
        print('Calculate rise velocity...')
        riseVelocity = get_riseVelocity(psi,vez,0.01)
        dat_riseVelocity.append(riseVelocity)
        print('done.')
        i = i + 1
    t_ = np.array(t_)
    dat_centerOfMass = np.array(dat_centerOfMass)
    dat_riseVelocity = np.array(dat_riseVelocity)
    dat_d_z = np.array(dat_d_z)
    dat_d_r = np.array(dat_d_r)
    return t_,dat_centerOfMass,dat_riseVelocity,dat_d_z,dat_d_r


def read_and_write_data(sub_dir,corr_list,start):
    # first check if the .dat file exist
    if os.path.exists('./' + sub_dir + '/'+sub_dir+'.dat'):
        print('The data file:  '+ './' + sub_dir + '/'+sub_dir+'.dat' + ' does exist already: reading data file ...')
        t_               = list()
        dat_centerOfMass = list()
        dat_riseVelocity = list()
        dat_d_z          = list()
        dat_d_r          = list()
        input = open('./' + sub_dir + '/'+sub_dir+'.dat', 'r')
        lines = input.readlines()
        for line in lines[1:]:
            numbers = line.split()
            t_.append(float(numbers[0]))
            dat_centerOfMass.append(float(numbers[1]))
            dat_riseVelocity.append(float(numbers[2]))
            dat_d_z.append(float(numbers[3]))
            dat_d_r.append(float(numbers[4]))
        t_ = np.array(t_)
        dat_centerOfMass = np.array(dat_centerOfMass)
        dat_riseVelocity = np.array(dat_riseVelocity)
        dat_d_z = np.array(dat_d_z)
        dat_d_r = np.array(dat_d_r)
    else:
        print('The data file:  '+ './' + sub_dir + '/'+sub_dir+'.dat' + ' does not exist yet: reading binary files ...')
        t_,dat_centerOfMass,dat_riseVelocity,dat_d_z,dat_d_r = read_bubble_data(sub_dir,corr_list,start)
        with open('./' + sub_dir + '/'+sub_dir+'.dat', 'w') as file:
            file.write('Time     CenterOfMass      BubbleVelocity       D_vertical     D_horizontal\n')
            for i in range(0,len(t_)):
                file.write(str(t_[i])+'     '+str(dat_centerOfMass[i])+'     '+str(dat_riseVelocity[i])+'     '+str(dat_d_z[i])+'     '+str(dat_d_r[i]))
                if i!=(len(t_)-1): file.write('\n')
    return t_,dat_centerOfMass,dat_riseVelocity,dat_d_z,dat_d_r

def FluTAS_to_VTK_2D(folder,scalarlist,vectorlist,z_corr,vel_corr,time):
    print('Setting Up VTK object for time: ' + str(time))
    iskip       = np.array([1,1,1]).astype(int)
    dummy,xp,yp,zp,xu,yv,zw = read_single_field_binary(folder+scalarlist[0],iskip,folder)
    zp = zp + z_corr
    z_grid,y_grid = np.meshgrid(zp,yp)
    dx = zp[1]-zp[0]
    imageVTK = vtk.vtkImageData()
    imageVTK.SetSpacing([dx,dx,dx])
    imageVTK.SetOrigin([0.0,yp[0],zp[0]])
    imageVTK.SetDimensions([1,len(yp),len(zp)])
    for filename in scalarlist:
        print('Add array '+filename+ ' to current VTK object')
        scalar,xp,yp,zp,xu,yv,zw = read_single_field_binary(folder+filename,iskip,folder)
        scalar = scalar[0]
        vtk_data = numpy_to_vtk(num_array=scalar.ravel(order = 'F'), deep=True, array_type=vtk.VTK_FLOAT)
        vtk_data.SetName(filename[0:3])
        imageVTK.GetPointData().AddArray(vtk_data)
    for vector in vectorlist:
        print('Add array vector fields:')
        print(vector[0])
        print(vector[1])
        print(vector[2])
        print('...to current VTK object')
        vex,xp,yp,zp,xu,yv,zw = read_single_field_binary(folder+vector[0],iskip,folder)
        vey,xp,yp,zp,xu,yv,zw = read_single_field_binary(folder+vector[1],iskip,folder)
        vez,xp,yp,zp,xu,yv,zw = read_single_field_binary(folder+vector[2],iskip,folder)
        print('Z velocity is corrected by ' + str(vel_corr))
        vez = vez - vel_corr
        vex = 0.0*vex[0]   # bounary condition, non penetration bc
        vey = vey[0]
        vez = vez[0]
        #assume homogenous neumann BC and interpolate vey and vez on scalar mesh!
        vey_int = vey
        vez_int = vez
        for j in range(1,len(vez)):
            for k in range(1,len(vez[0])):
                vez_int[j][k] = 0.5*(vez[j][k]+vez[j  ][k-1])
                vey_int[j][k] = 0.5*(vey[j][k]+vey[j-1][k  ])
        velocity = algs.make_vector(vex.ravel(order = 'F'),vey_int.ravel(order = 'F'),vez_int.ravel(order = 'F'))
        vtk_data = dsa.numpyTovtkDataArray(velocity)
        vtk_data.SetName(vector[3])
        imageVTK.GetPointData().AddArray(vtk_data)
    print('Write VTK object for time '+str(time))
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(folder+('000000'+'{:0.0f}'.format(1e6*time))[-7:-1]+'.vti')
    writer.SetInputData(imageVTK)
    writer.Update()
    writer.Write()
    return

def write_VTK_all(sub_dir,corr_list,start,selection):
    #create list of timestamps with times
    folder = './' + sub_dir + '/'
    files = glob.glob('./' + sub_dir + '/psi_fld_*.bin')
    files.sort()
    stamps = [[],[]]
    for file in files:
        time     = get_time(file,folder)
        filename = file.split('/')[-1].split('\\')[-1]
        stamp    = re.findall(r'\d+',file)[-1]
        stamps[0].append(stamp)
        stamps[1].append(time)
    for out_num in selection:
        stamp    = stamps[0][out_num]
        time     = stamps[1][out_num]
        z_corr   = get_corr_loc(time,corr_list)
        vel_corr = get_corr_vel(time,corr_list)
        scalarlist = ['pre_fld_'+stamp+'.bin','ent_fld_'+stamp+'.bin','s_r_fld_'+stamp+'.bin','psi_fld_'+stamp+'.bin']
        vectorlist = [['vex_fld_'+stamp+'.bin','vey_fld_'+stamp+'.bin','vez_fld_'+stamp+'.bin','vel']]
        FluTAS_to_VTK_2D(folder,scalarlist,vectorlist,z_corr,vel_corr,time+start)
    return

def postprocess_florschuetz(sub_dir,selection,start):
    corr_list = get_inflow_data('./' + sub_dir + '/')
    write_VTK_all(sub_dir,corr_list,start,selection)
    data = read_and_write_data(sub_dir,corr_list,start)
    return
    

################## problem specific data ##########
lamda_1 = 0.0199
lamda_2 = 0.1654
cp1 = 1806
cp2 = 3185
T_sat = 0.0
latent_heat = 849.9*1e3
T_wall = T_sat
T_far = 3.1
rho1 = 1.647
rho2 = 736.44
mu1 = 10.2*1e-6
mu2 = 448.7*1e-6
alpha1 = lamda_1/(rho1*cp1)
alpha2 = lamda_2/(rho2*cp2)
x_start = 0.00021
x0 = np.array([0])
res = minimize(opti, x0, method='nelder-mead',args=(rho1,rho2,cp1,cp2,T_far,T_sat,latent_heat),
               options={'xatol': 1e-20, 'disp': False})
eps_num = res.x[0]
start = (1.0/alpha2)*(x_start/(2.0*eps_num))**2
print('Epsilon value: '+str(eps_num))

################################################## 

# PostProcessing setting:
selection = [0,2,4]
#Postprocess
postprocess_florschuetz('data',selection,start)

