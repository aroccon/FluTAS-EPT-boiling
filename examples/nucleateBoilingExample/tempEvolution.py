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
params = {'legend.fontsize': 12,
          'legend.loc':'best',
          'figure.figsize': (15, 5),
         'lines.markerfacecolor':'none',
         'axes.labelsize': 12,
         'axes.titlesize': 12,
         'xtick.labelsize':12,
         'ytick.labelsize':12,
         'grid.alpha':0.6}
pylab.rcParams.update(params)




fig, ax = plt.subplots(figsize=(8, 2.5))

input = open('temps.txt', 'r')
outflow = list()

lines = input.readlines()
for line in lines[3:]:
    outflow.append(float(line.split()[-1]))
input.close()


ax.plot(outflow[:],c='hotpink')
ax.set_ylabel('Max. wall superheat [K]')
ax.set_xlabel('Timesteps')

fig.savefig('Evolution.pdf', format='pdf', bbox_inches='tight',transparent=True)
