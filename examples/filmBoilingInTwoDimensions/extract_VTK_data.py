import csv
import sys
import glob
import numpy as np
import copy
import os
import re
import math
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
        scalar = scalar.mean(axis=(0))
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
        vex = vex.mean(axis=(0))
        vey = vey.mean(axis=(0))
        vez = vez.mean(axis=(0))
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
    vtiname = ('000000'+'{:0.0f}'.format(1e3*time))[-7:]+'.vti'
    print('Printing file: '+ vtiname)
    writer.SetFileName(folder+vtiname)
    writer.SetInputData(imageVTK)
    writer.Update()
    writer.Write()
    return

def write_VTK_all(sub_dir,start,selection):
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
        z_corr   = 0.0
        vel_corr = 0.0
        scalarlist = ['pre_fld_'+stamp+'.bin','ent_fld_'+stamp+'.bin','s_r_fld_'+stamp+'.bin','psi_fld_'+stamp+'.bin']
        vectorlist = [['vex_fld_'+stamp+'.bin','vey_fld_'+stamp+'.bin','vez_fld_'+stamp+'.bin','vel']]
        FluTAS_to_VTK_2D(folder,scalarlist,vectorlist,z_corr,vel_corr,time+start)
    return

def postprocess_Sun2014_vti(sub_dir,selection,start):
    print('Write all VTK files...')
    write_VTK_all(sub_dir,start,selection)
    print('Writing VTK files: Done ...')
    return
    

# PostProcessing setting:
start = 0.0
selection = list(range(0,11))


#Extract the data

postprocess_Sun2014_vti('data',selection,start)
