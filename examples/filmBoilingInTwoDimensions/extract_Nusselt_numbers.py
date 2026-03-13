import csv
import sys
import glob
import numpy as np
import copy
import os
import re


#some subroutines

def read_and_write_film_data(sub_dir,start):
    # first check if the .dat file exist
    if os.path.exists('./' + sub_dir + '/data.dat'):
        print('The data file:  '+ './' + sub_dir + '/'+sub_dir+'.dat' + ' does exist already: reading data file ...')
        times            = list()
        Nusselts          = list()
        input = open('./' + sub_dir + '/'+sub_dir+'.dat', 'r')
        lines = input.readlines()
        for line in lines[1:]:
            numbers = line.split()
            times.append(float(numbers[0]))
            Nusselts.append(float(numbers[1]))
        times = np.array(times)
        Nusselts = np.array(Nusselts)
    else:
        print('The data file:  '+ './' + sub_dir + '/data.dat' + ' does not exist yet: reading binary files ...')
        times,Nusselts = get_Nusselt_vs_time(sub_dir)
        with open('./' + sub_dir + '/'+sub_dir+'.dat', 'w') as file:
            file.write('Time     Nusseltnumber\n')
            for i in range(0,len(times)):
                file.write(str(times[i])+'     '+str(Nusselts[i]))
                if i!=(len(times)-1): file.write('\n')
    return times,Nusselts

def get_Nusselt_vs_time(sub_dir):
    folder = './' + sub_dir + '/'
    times = list()
    Nusselts = list()
    with open(folder+'log.out', 'r') as fp:
        for l_no, line in enumerate(fp):
            # search string
            if (r'The Nusselt number is:' in line):
                time = float(line.split()[-1])
                Nusselt = float(line.split()[4])
                times.append(time)
                Nusselts.append(Nusselt)
    return times,Nusselts
    
def postprocess_Sun2014_Nuss(sub_dir,start):
    print('Extract Nusselt number vs time...')
    times,Nusselts = read_and_write_film_data(sub_dir,start)
    print('Extracting the nusseltnumber: Done..')
    return

# PostProcessing setting:
start = 0.0


#Extract the data

postprocess_Sun2014_Nuss('data',start)
