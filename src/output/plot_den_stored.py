#!/usr/bin python
import os
import matplotlib as mpl
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
try:
    import cPickle as cp
except:
    import pickle as cp


version = 'v0'
directory = 'plot_den_time_' + version + '/'

if not os.path.exists(directory):
    os.makedirs(directory)

dt = cp.load(open(version, 'r'))
total_xlat = []
for d in dt:
    if d.xlat not in total_xlat:
        total_xlat.append(d.xlat)
total_xlat.sort()
print(total_xlat)

i = 0
check = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
denAnalise = []

for d in dt:
    # plt.imshow(d.den, cmap=cm.rainbow, interpolation='nearest')
    if (d.xlat == 12):
        denAnalise.append(d)

for d in denAnalise:
    if i in check:
        plt.figure(i)
        plt.title()
        plt.imshow(d.den, cmap=cm.autumn)
        plt.colorbar()
        plt.savefig(directory + 'den-' + str(i) + '.png')

    i += 1
