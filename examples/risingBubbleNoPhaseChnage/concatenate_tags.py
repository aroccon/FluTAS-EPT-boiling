#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@(     #@@@@@@@@@@#        &@        @@                  @
#@@@@@@@@@@@@@@@@#/@@@@@      #@@@@@@@@         &@@@        @@                  @
#@@@@@@@@@@@@(//////@@@@@     #@@@@@@         @@@@@@        @@                  @
#@@@@@@@@@@@@@(//////@@@@#    #@@@&         @@@@@@@@        @@@@@@@        @@@@@@
#@@@@@@@/&@@@@@@#/////&@@@    #@,         @@@@@@@@@@        @@@@@@@        @@@@@@
#@@@@@(/////%@@@@@%////#@@@   #@        @@@@@@@@@@@@        @@@@@@@        @@@@@@
#@@@@//////////#@@@@&////@@/  #@         @@@@@@@@@@@        @@@@@@@        @@@@@@
#@@@@@@@@@////////(@@@&///@@  #@@#         @@@@@@@@@        @@@@@@@        @@@@@@
#@@@@@@@@@@@@@@%/////(@@@//@& #@@@@@         @@@@@@@        @@@@@@@        @@@@@@
#@&////////(#@@@@@@@@////@@/@ #@@@@@@@         &@@@@        @@@@@@@        @@@@@@
#@////////////////////(@@@#/@#%@@@@@@@@@&        &@@        @@@@@@@        @@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#       Written by Lorenz Weber (https://orcid.org/0000-0002-1657-2236)
#       in August 2025
#
#       This python script concatenates the fractional bubble information that
#       was distributed over many processes using MPI paralell computing. The
#       fractional bubble information is located in ./data/post/tagging. At the
#       end of script a list is dropped in ./data/post/concatenate containing
#       for each time step a list of all bubbles in the computational domain
#       including their bubble-averaged properties. This script produces a log
#       output for this process. If the bubble data exhibits a periodicity in 
#       the x-direction, this script will be able to account for that!
#       
#-------------------------------------------------------------------------------


import csv
import sys
import glob
import numpy as np
import copy
import math
import os
import re
from scipy import interpolate
from sympy.solvers import solve,nsolve
from scipy.special import erf,erfc
from sympy import Symbol,Eq
import matplotlib.cm as cm
import sympy
import scipy.integrate as integrate
from scipy.optimize import minimize
import struct
from math import floor

################################
#   required subroutines       #
################################

def read_binary_integers(filename):
    integers = []
    byte_order='@'
    format_char = 'i'
    byte_size = 4
    format_string = byte_order + format_char
    with open(filename, 'rb') as f:
        while(True):
            # Read the required number of bytes for one integer
            data = f.read(byte_size)
            if not data:  # End of file reached before reading all integers
                break
            # Unpack the bytes into an integer
            integer_value = struct.unpack(format_string, data)[0]
            integers.append(integer_value)
    return integers

def read_binary_integers6(filename,number_of_entries):
    integers = []
    byte_order='@'
    format_char = 'i'
    byte_size = 4
    format_string = byte_order + format_char
    with open(filename, 'rb') as f:
        while(True):
            # Read the required number of bytes for one integer
            data = f.read(byte_size)
            if not data:  # End of file reached before reading all integers
                break
            # Unpack the bytes into an integer
            integer_value = struct.unpack(format_string, data)[0]
            integers.append(integer_value)
    integers = list(np.reshape(integers, (number_of_entries,6)))
    return integers


def read_binary_floats(filename):
    floats = []
    with open(filename, 'rb') as f:
        while(True):
            # Read the required number of bytes for one integer
            data = f.read(8)
            if not data:  # End of file reached before reading all integers
                break
            # Unpack the bytes into an integer
            float_value = struct.unpack('d', data)[0]
            floats.append(float_value)
    return floats


def get_time(time_st,folder):
    stamp = ' ' + str(time_st) + ' '
    time = 0
    found = False
    if os.path.exists(folder+"log_visu_3d.out"):
        with open(folder+"log_visu_3d.out", 'r') as fp:
            for l_no, line in enumerate(fp):
                # search string
                if (int(line.split()[-1]) == time_st):
                    found = True
                    time = float(line.split()[-2])
                    break
    if not found:
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


################################
#   read in bubble tag info    #
################################

sub_dir = 'data/post/tagging'
bub_id_files = glob.glob('./' + sub_dir + '/dridfld*.bin')
bub_id_files.sort()

proc_id_files = glob.glob('./' + sub_dir + '/idmyfld*.bin')
proc_id_files.sort()

proc_id_share_files = glob.glob('./' + sub_dir + '/procfld*.bin')
proc_id_share_files.sort()

bub_id_share_files = glob.glob('./' + sub_dir + '/idshfld*.bin')
bub_id_share_files.sort()


################################
#     which files to load      #
################################


#                 w1      w2     w3     w4     w5
data_names    = ['vold','volx','voly','volz','vols','xpos','xmin','xmax','ypos','ymin','ymax','zpos','zmin','zmax','uvel','vvel','wvel','conl','drya','mtot','mmin','mmax','surf','pfmx']
data_action   = [   5  ,  5   ,  5   ,  5   ,   5  ,  9   ,  10  ,  11  ,  8   ,  6   ,  7   ,  8   ,  6   ,  7   ,  2   ,  3   ,  4   ,  5   ,  5   ,   5  ,   6  ,  7   ,   5  ,   7  ]
# Action 1 -> weighted sum with w1
# Action 2 -> weighted sum with w2
# Action 3 -> weighted sum with w3
# Action 4 -> weighted sum with w4
# Action 5 -> simple sum
# Action 6 -> min value
# Action 7 -> max value
# Action 8 -> weighted sum with w5
# Action 9 -> Action 8 with periodic correction
# Action 10-> Action 6 with periodic correction
# Action 11-> Action 7 with periodic correction

################################
#         load files.          #
################################

data_files = list()
for name in data_names:
    files = glob.glob('./' + sub_dir + '/'+name+'fld*.bin')
    files.sort()
    data_files.append(files)
data_files = [list(i) for i in zip(*data_files)] #transpose list

################################
#   input comp. settings       #
################################

# The following is important for accurately estimating the bubble radius from the bubble volume

c_angle = 179.0*np.pi/180.0

# The following is important for accounting for spanwise periodicity in x-direction for channelflows

procx = 4.0
procy = 8.0
xlen = 1 

################################
#   concatenate bubble info    #
################################


bubbes      = list()
for timestamp in range(0,len(bub_id_files)): # loop over all output times
    bubble_counter = 0
    time = int(bub_id_files[timestamp][-13:-4]) #get timestamp
    print('Processing timestep '+str(time)+' ('+str(timestamp+1)+'/'+str(len(bub_id_files))+')')
    time_r = get_time(time,'./data/')  # obtain real time from logfile
    # read binaries
    bub_ids        = read_binary_integers(bub_id_files[timestamp])
    proc_ids       = read_binary_integers(proc_id_files[timestamp])
    proc_id_shares = read_binary_integers6(proc_id_share_files[timestamp],len(bub_ids))
    bub_id_shares  = read_binary_integers6(bub_id_share_files[timestamp],len(bub_ids))
    data = list()
    for file in data_files[timestamp]:
        data.append(read_binary_floats(file))
    data = [list(i) for i in zip(*data)]
    #create a unique bubble fraction id
    indentifiers = [bub_ids,proc_ids]   # local bubble id plus process id is unique identifier
    indentifiers = [list(i) for i in zip(*indentifiers)] # tranpose data
    #loop until all bubble fractions have found their corresponding counterparts in other processors
    while len(bub_ids)!=0:
        bubble_counter   = bubble_counter + 1
        bubble_fractions = {(bub_ids[0],proc_ids[0])} # add the first bubble fraction to the bubble list
        added = True
        while(added):  # loop over all fractional bubbes infos and see if they belong one of the fractional infos that are already in the list
            added = False
            for fraction in bubble_fractions:
                fraction_list = list(fraction)
                frac_index = indentifiers.index(fraction_list)
                for direction in range(0,6):
                    share_id = (bub_id_shares[frac_index][direction],proc_id_shares[frac_index][direction])
                    if ((share_id[0]!=0)&(share_id[1]>=0)&(not(share_id in bubble_fractions))): #only add if the neighbor is a valid entry and bubble the bubble infor is not already in the list
                        added = True
                        break
                if(added):
                    break
            if(added):
                bubble_fractions.add(share_id)
        len_old = len(bubble_fractions)
        len_new = 0
        while(len_new!=len_old): # as long as not all corresponding tags are found
            len_new = len_old
            len_old = len(bubble_fractions)
            for entry in range(1,len(bub_ids)): #now add for all neighbor tags of all fractional bubble infos to the list
                #print(bubble_fractions,entry)
                frac_id = (bub_ids[entry],proc_ids[entry])
                for direction in range(0,6):
                    share_id = (bub_id_shares[entry][direction],proc_id_shares[entry][direction])
                    if (share_id in bubble_fractions):
                        bubble_fractions.add(frac_id)
                len_new = len(bubble_fractions)
                if (len_new!=len_old):  #if a new bubbe tag is added to the list, the list needs to be reevaluated
                    break
            added = False
            if (len_new!=len_old):
                added = True
            while(added): # loop over all fractional bubbes infos and see if they belong one of the fractional infos that are already in the list
                added = False
                for fraction in bubble_fractions:
                    fraction_list = list(fraction)
                    frac_index = indentifiers.index(fraction_list)
                    for direction in range(0,6):
                        share_id = (bub_id_shares[frac_index][direction],proc_id_shares[frac_index][direction])
                        if ((share_id[0]!=0)&(share_id[1]>=0)&(not(share_id in bubble_fractions))):
                            added = True
                            break
                    if(added):
                        break
                if(added):
                    print('Important')   #It is important that this tag was found: Therefore this is marked in the log file
                    bubble_fractions.add(share_id)
        # now all fractiopnal bubble info is found that belongs to that bubble
        bub_data = [time, time_r, bubble_counter]  # create general bubble info
        for data_id in range(0,len(data_names)):  
        #fill the bubble properties with generic start values:
            if(data_action[data_id]==6):
                bub_data.append(1.0e8)
            elif(data_action[data_id]==7):
                bub_data.append(-1.0e8)
            elif(data_action[data_id]==10):
                bub_data.append(1.0e8)
            elif(data_action[data_id]==11):
                bub_data.append(-1.0e8)
            else:
                bub_data.append(0.0)
        # now detect if the bubbes info is spread over a periodic direction
        sidelist = set()
        for fraction in bubble_fractions:
            proc = float(list(fraction)[1])
            fraction_list = list(fraction)
            frac_index = indentifiers.index(fraction_list)
            front_proc = proc_id_shares[frac_index][0]
            back_proc  = proc_id_shares[frac_index][1]
            if((floor(proc/procy)==0)&(floor(front_proc/procy)==(procx-1))):
                side = 1
            elif((floor(back_proc/procy)==0)&(floor(proc/procy)==(procx-1))):
                side = 1
            else:
                side = 0
            sidelist.add(side)
        periodic_correct = False
        # is correction required
        if(1 in sidelist):
            periodic_correct = True
        for fraction in bubble_fractions:
            proc = float(list(fraction)[1])
            fraction_list = list(fraction)
            frac_index = indentifiers.index(fraction_list)
            frac_data = data[frac_index]
            side_corr = 0.0
            if(periodic_correct&(floor(proc/procy)<floor(procx/2))):
                side_corr = xlen   #shift bubble info to one side
            for data_id in range(0,len(data_names)):  
            # combine bubble fractional info depending on data type
                match data_action[data_id]:
                    case 1:
                        bub_data[3+data_id] = bub_data[3+data_id] + frac_data[data_id]*frac_data[0]
                    case 2:
                        bub_data[3+data_id] = bub_data[3+data_id] + frac_data[data_id]*frac_data[1]
                    case 3:
                        bub_data[3+data_id] = bub_data[3+data_id] + frac_data[data_id]*frac_data[2]
                    case 4:
                        bub_data[3+data_id] = bub_data[3+data_id] + frac_data[data_id]*frac_data[3]
                    case 5:
                        bub_data[3+data_id] = bub_data[3+data_id] + frac_data[data_id]
                    case 6:
                        bub_data[3+data_id] = min(bub_data[3+data_id],frac_data[data_id])
                    case 7:
                        bub_data[3+data_id] = max(bub_data[3+data_id],frac_data[data_id])
                    case 8:
                        bub_data[3+data_id] = bub_data[3+data_id] + frac_data[data_id]*frac_data[4]
                    case 9:
                        bub_data[3+data_id] = bub_data[3+data_id] + (frac_data[data_id]+side_corr)*frac_data[4]
                    case 10:
                        bub_data[3+data_id] = min(bub_data[3+data_id],frac_data[data_id]+side_corr)
                    case 11:
                        bub_data[3+data_id] = max(bub_data[3+data_id],frac_data[data_id]+side_corr)
        if(bub_data[7]==0):# this is an error because it means that a fractional bubble info has not found its main bubble that it belongs to
            for fraction in bubble_fractions:
                # output info for degugging
                fraction_list = list(fraction)
                frac_index = indentifiers.index(fraction_list)
                print('Zero mass detected! At time: '+str(time_r))   
                print('Bubble count: '+str(bubble_counter)+'  '+str(periodic_correct))
                print(fraction_list,proc_id_shares[frac_index],bub_id_shares[frac_index])
                print(data[frac_index])
        # finish averaging by deviding by the sum of the weights
        for data_id in range(0,len(data_names)):
            match data_action[data_id]:
                case 1:
                    bub_data[3+data_id] = bub_data[3+data_id]/bub_data[3]
                case 2:
                    bub_data[3+data_id] = bub_data[3+data_id]/bub_data[4]
                case 3:
                    bub_data[3+data_id] = bub_data[3+data_id]/bub_data[5]
                case 4:
                    bub_data[3+data_id] = bub_data[3+data_id]/bub_data[6]
                case 8:
                    bub_data[3+data_id] = bub_data[3+data_id]/bub_data[7]
                case 9:
                    bub_data[3+data_id] = bub_data[3+data_id]/bub_data[7]
        for data_id in range(0,len(data_names)):
            pf_loc = 3+data_names.index('pfmx')
            # in case the maximum phasefile value inside a bubble does not reach 0.5, the empty data fields are filled with estimated info
            if(bub_data[pf_loc]<=0.5):
                match data_names[data_id]:
                    case 'xmin':
                        vol = bub_data[3+data_names.index('vols')]
                        rad = (3.0*vol/(4.0*np.pi))**(1.0/3.0)
                        if(bub_data[3+data_names.index('drya')]> 0.0):
                            rad1 = ((3.0*vol/np.pi)/((2+np.cos(c_angle))*((1-np.cos(c_angle))**2)))**(1.0/3.0)
                            area = bub_data[3+data_names.index('drya')]
                            rad2 = ((area*0.5/np.pi)/((1-np.cos(c_angle))))**(1.0/2.0)
                            rad  = 0.5*(rad1+rad2)
                        bub_data[3+data_id] = bub_data[3+data_names.index('xpos')]-rad
                    case 'xmax':
                        vol = bub_data[3+data_names.index('vols')]
                        rad = (3.0*vol/(4.0*np.pi))**(1.0/3.0)
                        if(bub_data[3+data_names.index('drya')]> 0.0):
                            rad1 = ((3.0*vol/np.pi)/((2+np.cos(c_angle))*((1-np.cos(c_angle))**2)))**(1.0/3.0)
                            area = bub_data[3+data_names.index('drya')]
                            rad2 = ((area*0.5/np.pi)/((1-np.cos(c_angle))))**(1.0/2.0)
                            rad  = 0.5*(rad1+rad2)
                        bub_data[3+data_id] = bub_data[3+data_names.index('xpos')]+rad
                    case 'ymin':
                        vol = bub_data[3+data_names.index('vols')]
                        rad = (3.0*vol/(4.0*np.pi))**(1.0/3.0)
                        if(bub_data[3+data_names.index('drya')]> 0.0):
                            rad1 = ((3.0*vol/np.pi)/((2+np.cos(c_angle))*((1-np.cos(c_angle))**2)))**(1.0/3.0)
                            area = bub_data[3+data_names.index('drya')]
                            rad2 = ((area*0.5/np.pi)/((1-np.cos(c_angle))))**(1.0/2.0)
                            rad  = 0.5*(rad1+rad2)
                        bub_data[3+data_id] = bub_data[3+data_names.index('ypos')]-rad
                    case 'ymax':
                        vol = bub_data[3+data_names.index('vols')]
                        rad = (3.0*vol/(4.0*np.pi))**(1.0/3.0)
                        if(bub_data[3+data_names.index('drya')]> 0.0):
                            rad1 = ((3.0*vol/np.pi)/((2+np.cos(c_angle))*((1-np.cos(c_angle))**2)))**(1.0/3.0)
                            area = bub_data[3+data_names.index('drya')]
                            rad2 = ((area*0.5/np.pi)/((1-np.cos(c_angle))))**(1.0/2.0)
                            rad  = 0.5*(rad1+rad2)
                        bub_data[3+data_id] = bub_data[3+data_names.index('ypos')]+rad
                    case 'zmin':
                        vol = bub_data[3+data_names.index('vols')]
                        rad = (3.0*vol/(4.0*np.pi))**(1.0/3.0)
                        if(bub_data[3+data_names.index('drya')]> 0.0):
                            rad1 = ((3.0*vol/np.pi)/((2+np.cos(c_angle))*((1-np.cos(c_angle))**2)))**(1.0/3.0)
                            area = bub_data[3+data_names.index('drya')]
                            rad2 = ((area*0.5/np.pi)/((1-np.cos(c_angle))))**(1.0/2.0)
                            rad  = 0.5*(rad1+rad2)
                        bub_data[3+data_id] = bub_data[3+data_names.index('zpos')]-rad
                    case 'zmax':
                        vol = bub_data[3+data_names.index('vols')]
                        rad = (3.0*vol/(4.0*np.pi))**(1.0/3.0)
                        if(bub_data[3+data_names.index('drya')]> 0.0):
                            rad1 = ((3.0*vol/np.pi)/((2+np.cos(c_angle))*((1-np.cos(c_angle))**2)))**(1.0/3.0)
                            area = bub_data[3+data_names.index('drya')]
                            rad2 = ((area*0.5/np.pi)/((1-np.cos(c_angle))))**(1.0/2.0)
                            rad  = 0.5*(rad1+rad2)
                        bub_data[3+data_id] = bub_data[3+data_names.index('zpos')]+rad
        # correct the x location of the bubble if the center of the bubble is outside the domain (for periodic correction)
        if((bub_data[3+data_names.index('xpos')]>xlen)&(periodic_correct)):
            print('X correction for bubble: '+str(bubble_counter))
            for data_id in range(0,len(data_names)):
                match data_names[data_id]:
                    case 'xmin':
                        bub_data[3+data_id] = bub_data[3+data_id]-xlen
                    case 'xmax':
                        bub_data[3+data_id] = bub_data[3+data_id]-xlen
                    case 'xpos':
                        bub_data[3+data_id] = bub_data[3+data_id]-xlen
        # finally append the data set for one bubble
        bubbes.append(bub_data)
        # and remove bubble fractional infos from the remaining list
        for fraction in bubble_fractions:
            fraction_list = list(fraction)
            frac_index = indentifiers.index(fraction_list)
            # create log info
            print('Bubble count: '+str(bubble_counter)+'  '+str(periodic_correct))
            print(fraction_list,proc_id_shares[frac_index],bub_id_shares[frac_index])
            print(data[frac_index])
            indentifiers.pop(frac_index)
            bub_ids.pop(frac_index)
            proc_ids.pop(frac_index)
            proc_id_shares.pop(frac_index)
            bub_id_shares.pop(frac_index)
            data.pop(frac_index)

################################
#   create output file         #
################################


with open('./data/post/concatenate/concatenated_bubbles.dat', 'w') as file:
    file.write('Timestamp  Time    bubble_ID')
    for name in data_names:
        file.write('      '+name)
    file.write('\n')
    for i in range(0,len(bubbes)):
        file.write(str(bubbes[i][0])+'     '+str(bubbes[i][1])+'     '+str(bubbes[i][2]))
        for data_id in range(0,len(data_names)):
            file.write('      '+str(bubbes[i][3+data_id]))
        if i!=(len(bubbes)-1): file.write('\n')


