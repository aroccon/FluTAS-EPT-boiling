# This is an example for how one can postprocess the avtivation frquency histogram of the nucleation sites

import csv
import sys
import glob
import numpy as np
import math
import os
import matplotlib.pyplot as plt
params = {'legend.fontsize': 12,
          'legend.loc':'best',
          'figure.figsize': (15, 5),
          'lines.markerfacecolor':'none',
          'axes.labelsize': 12,
          'axes.titlesize': 12,
          'xtick.labelsize':12,
          'ytick.labelsize':12}
import matplotlib.pylab as pylab
pylab.rcParams.update(params)



files = ['log.out']
input = open('randomNucleationSites.in', 'r')
lines = input.readlines()
sites = list()
tag = 0
for line in lines[2:]:
    tag = tag + 1
    numbers = line.split()
    x = float(numbers[0])
    y = float(numbers[1])
    z = float(numbers[2])
    T = float(numbers[3])
    sites.append([tag,x,y,z,T])


events = list()
for _ in range(0,len(sites)):
    events.append(list())
    
min_time = 1e10
max_time = 0.0

events_count = 0

for file in files:
    input = open(file, 'r')
    lines = input.readlines()
    for line in lines:
        if 'has been activated at surface temp:' in line:
            events_count = events_count + 1
            rest = line.strip().split('Nucleation site no.')[1].strip()
            nuc_id = int(rest.split()[0])-1
            rest = rest.split('has been activated at surface temp:')[1].strip()
            temp = float(rest.split()[0])
            rest = rest.split('at time')[-1].strip()
            time = float(rest)
            min_time = min(min_time,time)
            max_time = max(max_time,time)
            events[nuc_id].append(time)
    input.close()



freqency = np.zeros(len(sites))
freq_no_zero = list()
zeros = 0
for i in range(0,len(sites)):
    freqency[i] = float(len(events[i]))/(max_time-min_time)
    if(len(events[i])==0):
        zeros = zeros + 1
    else:
        freq_no_zero.append(float(len(events[i]))/(max_time-min_time))

freq_no_zero = np.array(freq_no_zero)


#Mean
mean_Weber = events_count/((max_time-min_time)*(len(sites)-zeros))  #

#36 bins
#
bins = np.logspace(np.log10(2),np.log10(10000),num = 37,base=10)
bins_ext = np.zeros(len(bins)+1)
bins_ext[1:] = bins

hist_weber, _ = np.histogram(freq_no_zero,bins = bins_ext,density = True,)
hist_weber

# Comparison data from A. Kossolapov  1 MPa, 500 kg/m²s, 1 MW/m² // Thesis p. 173,  Fig 5-3
frequencies_Koss = [2.1975575789572264, 2.716366491606049, 3.5752793253452575, 4.56031492186923, 5.7261314584045, 7.303754576385493, 9.316033224102007, 11.697619972746635, 14.80379888473144, 18.442953311699963, 23.896467278263955, 30.005452555612948, 38.27234198368184, 48.43515645867369, 61.29659852476358, 78.18460250943131, 98.94569443432881, 124.24055427774162, 159.71931041387796, 198.98242389514218, 251.8200960368871, 326.282379142459, 397.030018497482, 498.52830732082145, 645.9413078719249, 811.0722409873163, 1034.5331110405464, 1299.0051549909642, 1618.3340243743487, 2064.2059342343396, 2716.906140142946, 3411.4665292080417, 4385.662030507793, 5421.048120538096, 7024.032259744001, 8614.407461887622]
weights_Koss = [0, 0, 0.6816578392288227, 0, 0, 0.1073418882701355, 0, 0.02265920952388175, 0.012328467394420685, 0.004676598326676731, 0.0024877354553875233, 0.0026617963498590845, 0.0011301780105701623, 0.0007703947006291849, 0.0006882753845803487, 0.00022297272846923507, 0.00034999024845077737, 0.00010596802134320273, 0.00010129649200722846, 0.00005388507993839957, 0.000038425133462543194, 0.00002560894034079315, 0.00001457593300969297, 0.000009286074590769966, 0.000007580891478870158, 0.000004513832659768975, 0.000011634139467879642, 0.000006329925540461916, 0.0000021452036709325255, 0, 0, 0,0, 0, 0, 0]
mean_Koss = 20.9

A, B = np.histogram(frequencies_Koss, bins=bins, weights=weights_Koss,density=False,)


fig, ax = plt.subplots(figsize=(5, 5))


ax.hist(frequencies_Koss, bins, weights=A, histtype='step', edgecolor='b',label = 'A. Kossolapov, PhD Thesis (2021)', fill=False)
ax.hist(bins_ext[:-1], bins_ext, histtype='step', edgecolor='r',label = 'Present (Simulation)', fill=False, weights=hist_weber)

ax.plot([mean_Koss,mean_Koss],[0,10],linestyle='dashed',c='b')
ax.plot([mean_Weber,mean_Weber],[0,10],linestyle='dashed',c='r')

lgd = ax.legend()
frame = lgd.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
ax.set_yscale('log', base=10)
ax.set_xscale('log', base=10)
ax.set_xlim([2,10000])
ax.set_ylim([0.0000001,10])
#ax.ticklabel_format(axis = 'y',style = 'sci',useMathText=True, scilimits=(0,0))
ax.set_ylabel(r'PDF')
ax.set_xlabel(r'Bubble departure frequency [Hz]')

fig.savefig('BubbleDepartureFrequencyPDF.pdf', format='pdf', bbox_inches='tight',transparent=True)


