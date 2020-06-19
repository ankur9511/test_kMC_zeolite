#!/usr/bin/env python2
# coding: utf-8

# In[ ]:


import numpy as np
import copy
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from matplotlib.animation import FuncAnimation, writers
import mpl_toolkits.mplot3d.axes3d as p3
import gc
from multiprocessing import Pool
import glob
import os
import sys
import time
from functools import partial
import copy


# In[ ]:


def writefile(t,x,y,z,fname,d=None):
	trajfile = open(fname,'ab')
	trajfile.write('%.8E' %t)
	trajfile.write('\t')
	trajfile.write('%.5f' %x)
	trajfile.write('\t')
	trajfile.write('%.5f' %y)
	trajfile.write('\t')
	trajfile.write('%.5f' %z)
	if (d!=None):
		trajfile.write('\t')
		trajfile.write('%.5f' %d)
	trajfile.write('\n')
	trajfile.close()


# In[ ]:


global ID,ndiff,nlist,d_alist,label,rate_list,sys_length,universal_site,site_nsim,dtraj_0_nsim,dtraj_k_nsim,tsteps,exitcondition,itr,start,minx,maxx,miny,maxy,minz,maxz


# In[ ]:


def kmc(): #ndiff,nlist,alist,d_alist,label,rate_list,sys_length,universal_site,site_nsim,dtraj_0_nsim,dtraj_k_nsim,tsteps,exitcondition,itr,start):
	#### Do diffusive kmc for every run
	#### dtraj: Trajectory
	#### deltraj: instantaneous displacement
	#### dt: time evolution
	#### pos: get initial position of particles
	#### How to assign initial positions of diffusing particles
	#### 1. Pick molecule that would diffuse at a time step
	#### 2. Find the last position of that molecule to use as new initial position
	#### 3. Reinitialize rate catalog to empty
	#### 4. Find the unoccupied neighbors
	#### 5. Pick the neighbor event that would take place
	global ID,ndiff,nlist,d_alist,label,rate_list,sys_length,universal_site,site_nsim,dtraj_0_nsim,dtraj_k_nsim,tsteps,exitcondition,itr,start
	time.sleep(int(itr)%10+1)
	itr = int(itr)
	run = str(itr)
	dtraj_0 = copy.deepcopy(dtraj_0_nsim[int(itr)])#np.zeros((ndiff,3))
	dtraj_k = copy.deepcopy(dtraj_k_nsim[int(itr)])#np.zeros((ndiff,3))
	site = copy.deepcopy(site_nsim[int(itr)])
	dtraj_kp1 = np.zeros((ndiff,3))
	dist = np.zeros(ndiff)
	dt = 0
	dt_0 = 0
	dt_f = 0
	t = 1
	universal_site+=[site]
	flag = True
	for p in glob.glob("traj"+run+'_*'):
		if os.path.exists(p):
			os.remove(p)
	for p in glob.glob("eo_traj_"+run):
		if os.path.exists(p):
			os.remove(p)
	for event_s in range(len(site)):
		dt = 0.0
		trajfile = open('traj'+run+'_'+str(event_s),'ab')
		trajfile.write('%.8E' %dt)
		trajfile.write('\t')
		trajfile.write('%.5f' %dtraj_0[event_s,0])
		trajfile.write('\t')
		trajfile.write('%.5f' %dtraj_0[event_s,1])
		trajfile.write('\t')
		trajfile.write('%.5f' %dtraj_0[event_s,2])
		trajfile.write('\t')
		trajfile.write('%.5f' %dist[event_s])
		trajfile.write('\n')
		trajfile.close()
	writeflag = -1
	del event_s
	#ID = {}
	#i = 0
	#while i < nlist.shape[0]:
	#	ID[nlist[i,0]] = i
	#	i = i+1
	#
	while flag:
		if (exitcondition == "timestep"):
			if (t > tsteps):
				flag = False
				f = open("eo_traj_"+run,'ab')
				f.write("tstep exit\n")
				f.close()
		elif (exitcondition == "surface_exit"):
			for j in site:
				if (label[j-1] == "ZGI" or label[j-1] == "GI"):
					flag = False
					f = open("eo_traj_"+run,'ab')
					f.write('nsteps = '+str(t)+'\n')
					f.write("surface exit\n")
					f.close()
		elif (exitcondition == "smooth boundary"):
			for j in site:
				if (dtraj_k[0,0] < minx or dtraj_k[0,0] > maxx or dtraj_k[0,1] < miny or dtraj_k[0,1] > maxy or dtraj_k[0,2] < minz or dtraj_k[0,2] > maxz):
					flag = False
					f = open("eo_traj_"+run,'ab')
					f.write('nsteps = '+str(t)+'\n')
					f.write("boundary exit\n")
					f.close()

		if (flag == False):
			if (writeflag==0):
				writefile(dt,dtraj_k[event_s,0],dtraj_k[event_s,1],dtraj_k[event_s,2],'traj'+run+'_'+str(event_s),dist[event_s])
			f = open('eo_traj_'+run,'ab')
			f.write("Ending simulation\n")
			f.close()
			break
		dtraj_kp1[:,:] = 0.0
		RatCat = []					#.... 3
		tot_rate = 0					#.... 3
		event_site = []
		event_nsite = []
		event_asite = []
		for j in range(len(site)):
			i = ID[site[j]]
			for k in range(len(nlist[i,1:])-len(nlist[i,nlist[i]==0])):
				if (writeflag == 0 and label[int(nlist[i,1+k])-1]=='ZGI'):
					writeflag = 1
					writefile(dt,dtraj_k[event_s,0],dtraj_k[event_s,1],dtraj_k[event_s,2],'traj'+run+'_'+str(event_s),dist[event_s])
				if (nlist[i,1+k] not in site):									#.... 4
					if (label[int(nlist[i,1+k])-1] in rate_list[label[int(nlist[i,0])-1]]):
						RatCat += [rate_list[label[int(nlist[i,0])-1]][label[int(nlist[i,1+k])-1]]]
						tot_rate += float(rate_list[label[int(nlist[i,0])-1]][label[int(nlist[i,1+k])-1]])
					else:
						RatCat+=[0]
				else:
					RatCat+=[0]
				event_site += [j]
				event_asite += [i]
				event_nsite += [k]
		pick_event = (1.0-np.random.random())*tot_rate
		for i in range(len(RatCat)):																#.... 5
			if pick_event < sum(RatCat[:i+1]) :														#.... 12
				event_s = event_site[i]
				event_a = event_asite[i]
				event_n = event_nsite[i]
				break
		if (exitcondition == "smooth boundary"):
			if (t==1 and (label[nlist[event_a,1+event_n]-1]=='ZGI' or d_alist[0,int(nlist[event_a,1+event_n])-1] < minx or d_alist[0,int(nlist[event_a,1+event_n])-1] > maxx or d_alist[1,int(nlist[event_a,1+event_n])-1] < miny or d_alist[1,int(nlist[event_a,1+event_n])-1] > maxy or d_alist[2,int(nlist[event_a,1+event_n])-1] < minz or d_alist[2,int(nlist[event_a,1+event_n])-1] > maxz)):
				continue
		elif (exitcondition == "surface_exit"):
			if (t==1 and label[nlist[event_a,1+event_n]-1]=='ZGI'):
				continue

		##### Update dt, dtraj_kp1 #####
		dt = dt-1.0*np.log(1.0-np.random.random())/tot_rate												#.... 14
		#print dt,1.0/tot_rate,RatCat 
		dtraj_kp1[event_s,:] =  d_alist[:,int(nlist[event_a,1+event_n])-1] - d_alist[:,int(nlist[event_a,0])-1]		#.... 15
		delta = d_alist[:,int(nlist[event_a,1+event_n])-1] - d_alist[:,int(nlist[event_a,0])-1]
		Lx = sys_length[label[nlist[event_a,0]-1]][label[nlist[event_a,1+event_n]-1]]['x']
		Ly = sys_length[label[nlist[event_a,0]-1]][label[nlist[event_a,1+event_n]-1]]['y']
		Lz = sys_length[label[nlist[event_a,0]-1]][label[nlist[event_a,1+event_n]-1]]['z']
		dtraj_kp1[event_s,0] = dtraj_kp1[event_s,0] - Lx*np.around(dtraj_kp1[event_s,0]/Lx)	#.... 16
		dtraj_kp1[event_s,1] = dtraj_kp1[event_s,1] - Ly*np.around(dtraj_kp1[event_s,1]/Ly)	
		dtraj_kp1[event_s,2] = dtraj_kp1[event_s,2] - Lz*np.around(dtraj_kp1[event_s,2]/Lz)
		dist[event_s] = dist[event_s]+np.sqrt(np.dot([dtraj_kp1[event_s,0],dtraj_kp1[event_s,1],dtraj_kp1[event_s,2]],[dtraj_kp1[event_s,0],dtraj_kp1[event_s,1],dtraj_kp1[event_s,2]]))
		dtraj_kp1 = dtraj_kp1 + dtraj_k	
		##### Updated dt, dtraj_kp1
		##### Update site, dtraj_k,t ######																
		dtraj_k[:,:] = dtraj_kp1[:,:]																	
		site[event_s] = nlist[event_a,1+event_n]																		
		t = t+1
		##### Updated site, dtraj_k, t ###
		writeflag = 0
		if (t%1000 == 0):
			writeflag = 1
			writefile(dt,dtraj_k[event_s,0],dtraj_k[event_s,1],dtraj_k[event_s,2],'traj'+run+'_'+str(event_s),dist[event_s])
		elif (t in [1,2]):
			writeflag = 1
			writefile(dt,dtraj_k[event_s,0],dtraj_k[event_s,1],dtraj_k[event_s,2],'traj'+run+'_'+str(event_s),dist[event_s])
	for event_s in range(len(site)):
		trajfile = open('traj'+run+'_'+str(event_s),'ab')
		trajfile.write('Finished')
	return universal_site


# In[ ]:


def animate(i,simstep,dlist,points):
	for j,k in zip(simstep[:,i],points):
		data = dlist[:,j-1:j]
		k.set_data(data[0:2])
		k.set_3d_properties(data[2])
	return points


# In[ ]:


def writepdb(fname,Lx,Ly,Lz,neighbors,dist,label):
	f = open(fname,'w')
	i = 1
	newID = {}
	for ID in neighbors:
		newID[ID[0][0]] = i
		i = i+1
	f.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%-11s%4d\n" % (10*Lx,10*Ly,10*Lz,90,90,90," P 1",1))
	
	i = 0
	for ID in neighbors:
		f.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format("ATOM",newID[ID[0][0]],label[ID[0][0]-1],"X","Res","X",1,"X",(10*dist[0,ID[0][0]-1]),(10*dist[1,ID[0][0]-1]),(10*dist[2,ID[0][0]-1]),1,0,"X","X"))
		i = i+1
	
	i = 0
	for ID in neighbors:
		f.write("CONECT%5d" % (newID[ID[0][0]]))
		for j in ID[1]:
			f.write("%5d" % (newID[j]))
		f.write("\n")
		i = i+1
	
	f.write("END")
	f.close()
	return newID


# In[ ]:


def show(siteID,neighbors):
	siteID,neighbors = mergelist(siteID,neighbors)
	print "Length of site list ",len(siteID)
	print "Length of neighbors ", len(neighbors)
	l = []
	ID = []
	aID = []
	for i in neighbors:
		for j in i:
			l+=[len(j)]
			for k in j:
				ID+=[k]
	for i in siteID:
		for j in i:
			aID += [j]
	print "Neighbors ", np.unique(np.asarray(l)), "Maximum site number ", max(ID), "Number of unique",len(np.unique(ID)) ,"Maximum site number ", max(aID)
	del ID,aID


# In[ ]:


def find_smooth_boundary():
	global nlist,label,d_alist,minx,maxx,miny,maxy,minz,maxz
	fin_size = 30
	size_wo_fins = 300
	cell_size_x = 2.01
	cell_size_y = 1.99
	cell_size_z = 1.34
	system_size_x = size_wo_fins+2*fin_size+6*cell_size_x # float(sys.argv[4])
	system_size_y = size_wo_fins+2*fin_size+6*cell_size_y # float(sys.argv[5])
	system_size_z = size_wo_fins+2*fin_size+6*cell_size_z # float(sys.argv[6])
	nx = int(system_size_x/cell_size_x)+1
	ny = int(system_size_y/cell_size_y)+1
	nz = int(system_size_z/cell_size_z)+1
	lx = int((nx-int(size_wo_fins/cell_size_x)+1)*0.5)
	ly = int((ny-int(size_wo_fins/cell_size_y)+1)*0.5)
	lz = int((nz-int(size_wo_fins/cell_size_z)+1)*0.5)
	#Smooth part in rough zeolite created according to: system[lx:-lx,ly:-ly,lz:-lz] = 1
	minx = cell_size_x * lx #Index i has left edge at i, right edge at i+1
	maxx = cell_size_x * (nx-lx+1)
	miny = cell_size_y * ly #Index i has left edge at i, right edge at i+1
	maxy = cell_size_y * (ny-ly+1)
	minz = cell_size_z * lz #Index i has left edge at i, right edge at i+1
	maxz = cell_size_z * (nz-lz+1)


# In[ ]:


def main(i):
	La = 2.01
	Lb = 1.99
	Lc = 1.34
	global	ID,ndiff,nlist,d_alist,label,rate_list,sys_length,universal_site,site_nsim,dtraj_0_nsim,dtraj_k_nsim,tsteps,exitcondition,itr,start,minx,maxx,miny,maxy,minz,maxz
	#alist = np.load(sys.argv[1]+"_sitelist.npy",allow_pickle=True)
	shape = np.load(sys.argv[1]+'_shape.npy',allow_pickle=True)
	print "Neighbor reading started"
	nlist = np.memmap(sys.argv[1]+"_neighbors.npy",dtype=int,shape=(shape[0],shape[1]),mode='r')
	print "Neighbor reading ended"
	label = np.load(sys.argv[1]+"_label.npy",allow_pickle=True,mmap_mode='r')
	print "Label read"
	if (sys.argv[4] == 'HighDesorption'):
		rate_list = np.load("FSHdes_ratelist.npy",allow_pickle=True)
	elif (sys.argv[4] == 'LowDesorption'):
		rate_list = np.load("FSLdes_ratelist.npy",allow_pickle=True)
	rate_list = rate_list.item()
	print "Ratelist read"
	sys_length = np.load(sys.argv[1]+"_length.npy",allow_pickle=True)
	sys_length = sys_length.item()
	print "PBCLength read"
	ndiff = 1 #number of guest molecules
	tsteps = 10000 # 0: No tstep condition, > 0 : apply tstep condition
	d_alist = np.load(sys.argv[1]+"_coordinates.npy",allow_pickle=True,mmap_mode='r')
	print "Coordinates read"
	type_exits = ["timestep","surface_exit","smooth boundary"]
	exitcondition = type_exits[1]
	find_smooth_boundary()
	universal_site=[]
	site_nsim = np.load('site_nsim.npy',allow_pickle=True,mmap_mode='r')
	dtraj_0_nsim =  np.load('dtraj_0_nsim.npy',allow_pickle=True,mmap_mode='r')
	dtraj_k_nsim = np.load('dtraj_k_nsim.npy',allow_pickle=True,mmap_mode='r')
	print "Initial positions read and KMC started"
	
	start = int(sys.argv[2])
	last = int(sys.argv[3])
	ID = {}
	i = 0
	while i < nlist.shape[0]:
		ID[nlist[i,0]] = i
		i = i+1
	print "Initial positions read and KMC started"
	f = open(sys.argv[5],'w',0)
	f.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%-11s%4d\n" % (365,365,365,90,90,90," P 1",1))
	for itr in range(start,last,1):
		f.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format("ATOM",itr,'Si',"X","Res","X",1,"X",(dtraj_0_nsim[int(itr),0,0]),(dtraj_0_nsim[int(itr),0,1]),(dtraj_0_nsim[int(itr),0,2]),1,0,"X","X"))
		kmc()#ndiff,neighbors,siteID,dlist,label,rate_list,wlength,universal_site,site,dtraj_0,dtraj_k,tsteps,exit_condition,j,start)
	f.close()
main(1)

