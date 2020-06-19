#!/usr/bin/env python2
# coding: utf-8

import numpy as np
import sys
from matplotlib import rcParams
params = {'backend': 'agg'}
rcParams['font.family'] = 'serif'
rcParams.update(params)
from matplotlib import pyplot as plt

##### How to use #####
##### python smooth.py size_wo_fins(nm) cell_x(nm) cell_y(nm) cell_z(nm) name_of_file
##### 

size_wo_fins = float(sys.argv[1])
cell_size_x = 2.01 #float(sys.argv[2])
cell_size_y = 1.99 #float(sys.argv[3])
cell_size_z = 1.34 #float(sys.argv[4])
name = sys.argv[2]

print "Smooth Cells in x: ",int(size_wo_fins/cell_size_x)
print "Smooth Cells in y: ",int(size_wo_fins/cell_size_y)
print "Smooth Cells in z: ",int(size_wo_fins/cell_size_z)

system = np.zeros((int(size_wo_fins/cell_size_x)+4,int(size_wo_fins/cell_size_y)+4,int(size_wo_fins/cell_size_z)+4),dtype=int)
system[2:-2,2:-2,2:-2] = 1
np.save(name+'_z',system)

ig_system = np.zeros(system.shape,dtype=int)
ig_system[1:-1,1:-1,1:-1] = 1
ig_system = ig_system-system
np.save(name+'_ig',ig_system)

g_system = np.ones(system.shape,dtype=int)
g_system = g_system - ig_system - system
np.save(name+'_g',g_system)

boundary = 0*g_system
boundary[0,:,:] = 1
boundary[:,0,:] = 1
boundary[:,:,0] = 1
boundary[-1,:,:]= 1
boundary[:,-1,:]= 1
boundary[:,:,-1]= 1
np.save(name+"_boundary",boundary)


xy = np.average(system,axis=2)
xz = np.average(system,axis=1)
yz = np.average(system,axis=0)

plt.imshow(xy,extent=[0,cell_size_x*np.size(system,axis=0),0,cell_size_y*np.size(system,axis=1)],aspect='auto')
plt.savefig(name+"xy")
plt.imshow(yz,extent=[0,cell_size_y*np.size(system,axis=1),0,cell_size_z*np.size(system,axis=2)],aspect='auto')
plt.savefig(name+"yz")
plt.imshow(xz,extent=[0,cell_size_x*np.size(system,axis=0),0,cell_size_z*np.size(system,axis=2)],aspect='auto')
plt.savefig(name+"xz")

