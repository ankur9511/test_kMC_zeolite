#!/usr/bin/env python2
# coding: utf-8

# In[1]:


import numpy as np
import sys
from matplotlib import rcParams
params = {'backend': 'agg'}
rcParams['font.family'] = 'serif'
rcParams.update(params)
from matplotlib import pyplot as plt


# <img src="kMC_finmatillustrate.png" />

# ### Input parameters
# ##### surf_fin_size: L2 (nm)
# ##### h_fin_size: L1 (nm)
# ##### d_fin_size: L3 (nm)
# ##### size_wo_fins: Core of the finned zeolite
# ##### cell_size_x: unit cell size in x (nm)
# ##### cell_size_y: unit cell size in y (nm)
# ##### cell_size_z: unit cell size in z (nm)
# ##### name: prefix name of files generated
# In[5]:


surf_fin_size = float(sys.argv[1])
h_fin_size = float(sys.argv[2])
d_fin_size = float(sys.argv[3])
size_wo_fins = float(sys.argv[4])
cell_size_x = 2.01 #float(sys.argv[5])
cell_size_y = 1.99 #float(sys.argv[6])
cell_size_z = 1.34 #float(sys.argv[7])
name = sys.argv[5]

# ### Total cuboidal system size = zeolite core + 2 * height of fins + 2 * three layers of unit cells

# In[ ]:


system_size_x = size_wo_fins+2*h_fin_size+6*cell_size_x # float(sys.argv[4])
system_size_y = size_wo_fins+2*h_fin_size+6*cell_size_y # float(sys.argv[5])
system_size_z = size_wo_fins+2*h_fin_size+6*cell_size_z # float(sys.argv[6])


# ### System size in terms of number of unit cells in x, y, and z 

# In[ ]:


system = np.zeros((int(system_size_x/cell_size_x)+1,int(system_size_y/cell_size_y)+1,int(system_size_z/cell_size_z)+1),dtype=int)
ig_system = np.zeros((int(system_size_x/cell_size_x)+1,int(system_size_y/cell_size_y)+1,int(system_size_z/cell_size_z)+1),dtype=int)
g_system = np.ones((int(system_size_x/cell_size_x)+1,int(system_size_y/cell_size_y)+1,int(system_size_z/cell_size_z)+1),dtype=int)


# ##### l(x,y,z): Half of the number of unit cells that make the zeolite core
# ##### system: Matrix with value 1 where unit cell is a zeolite
# ##### ig_system: matrix with value 1 where unit cells comprise of the layer between zeolite and gas (not of use in current paper description
# ##### g_system: matrix with value 1 where unit cells comprise the external gaseous region
# #### Converting various lengths into analogous number of unit cells

# In[ ]:


lx = int((np.size(system,axis=0)-int(size_wo_fins/cell_size_x)+1)*0.5)
ly = int((np.size(system,axis=1)-int(size_wo_fins/cell_size_y)+1)*0.5)
lz = int((np.size(system,axis=2)-int(size_wo_fins/cell_size_z)+1)*0.5)
system[lx:-lx,ly:-ly,lz:-lz] = 1
ig_system[lx-1:-lx+1,ly-1:-ly+1,lz-1:-lz+1] = 1
sfx = int(surf_fin_size/cell_size_x)
sfy = int(surf_fin_size/cell_size_y)
sfz = int(surf_fin_size/cell_size_z)
hfx = int(h_fin_size/cell_size_x)
hfy = int(h_fin_size/cell_size_y)
hfz = int(h_fin_size/cell_size_z)


# In[ ]:


dx = int(d_fin_size/cell_size_x)
dy = int(d_fin_size/cell_size_y)
dz = int(d_fin_size/cell_size_z)


# In[ ]:


fx = int((np.size(system,axis=0)-2*lx)/(sfx+dx))
fy = int((np.size(system,axis=1)-2*ly)/(sfy+dy))
fz = int((np.size(system,axis=2)-2*lz)/(sfz+dz))


# In[ ]:


stx = int(lx+0.5*((np.size(system,axis=0)-2*lx)-fx*sfx-(fx-1)*dx))
lax = int(stx+(fx-1)*sfx+(fx-1)*dx)


# In[ ]:


sty = int(ly+0.5*((np.size(system,axis=1)-2*ly)-fy*sfy-(fy-1)*dy))
lay = int(sty+(fy-1)*sfy+(fy-1)*dy)


# In[ ]:


stz = int(lz+0.5*((np.size(system,axis=2)-2*lz)-fz*sfz-(fz-1)*dz))
laz = int(stz+(fz-1)*sfz+(fz-1)*dz)


# ### Decorate fins with a cubic arrangment on the zeolite core in a symetrical fashion

# In[ ]:


for i in range(sty,lay+sfy+dy,sfy+dy):
	for j in range(stz,laz+sfz+dz,sfz+dz):
		system[lx-hfx:-lx+hfx,i:i+sfy,j:j+sfz] = 1
		ig_system[lx-hfx-1:-lx+hfx+1,i-1:i+sfy+1,j-1:j+sfz+1] = 1


# In[ ]:


for i in range(stx,lax+sfx+dx,sfx+dx):
	for j in range(stz,laz+sfz+dz,sfz+dz):
		system[i:i+sfx,ly-hfy:-ly+hfy,j:j+sfz] = 1
		ig_system[i-1:i+sfx+1,ly-hfy-1:-ly+hfy+1,j-1:j+sfz+1] = 1


# In[6]:


for i in range(sty,lay+sfy+dy,sfy+dy):
	for j in range(stx,lax+sfx+dx,sfx+dx):
		system[j:j+sfx,i:i+sfy,lz-hfz:-lz+hfz] = 1
		ig_system[j-1:j+sfx+1,i-1:i+sfy+1,lz-hfz-1:-lz+hfz+1] = 1
ig_system = ig_system - system
g_system = g_system - ig_system - system

np.save(name+"_z",system)
np.save(name+"_ig",ig_system)
np.save(name+"_g",g_system)
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


# In[ ]:


plt.imshow(xy)
plt.savefig(name+"xy_z")
plt.imshow(yz)
plt.savefig(name+"yz_z")
plt.imshow(xz)
plt.savefig(name+"xz_z")


# In[ ]:


xy = np.average(ig_system,axis=2)
xz = np.average(ig_system,axis=1)
yz = np.average(ig_system,axis=0)


# In[ ]:


plt.imshow(xy)
plt.savefig(name+"xy_ig")
plt.imshow(yz)
plt.savefig(name+"yz_ig")
plt.imshow(xz)
plt.savefig(name+"xz_ig")


# In[ ]:


xy = np.average(g_system,axis=2)
xz = np.average(g_system,axis=1)
yz = np.average(g_system,axis=0)


# In[ ]:


plt.imshow(xy)
plt.savefig(name+"xy_g")
plt.imshow(yz)
plt.savefig(name+"yz_g")
plt.imshow(xz)
plt.savefig(name+"xz_g")

