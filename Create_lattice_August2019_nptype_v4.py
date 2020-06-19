#!/usr/bin/env python2
# coding: utf-8
#
# Scroll to bottom to understand each component
#
#
# In[1]:


import numpy as np
import copy
import glob
import os
import sys
import time
from multiprocessing import Pool
import gc
import datetime


# In[2]:


global siteID,neighbors,label,Zeo,Zeo_num,ZG,ZG_num,Gas,Gas_num,nZuc,nZGuc,nGuc,num_neighbors
t0 =  time.time()


# In[3]:


def readlattice(lfile,numstart=0):
	xyz = np.load(lfile)
	xyz_num = np.load(lfile)
	k = 0
	print "Lattice reading started"
	for i in np.ndindex(xyz.shape):
		if (xyz[i] == 1):
			xyz_num[i] = int(k+1+numstart)
			k = k+1
	return xyz,xyz_num,k


# In[6]:


def initialize_z(lfile_z,lfile_ig,lfile_g):
	global Zeo,Zeo_num,nZuc,ZG,ZG_num,nZGuc,Gas,Gas_num,nGuc,siteID,neighbors,num_neighbors
	Zeo,Zeo_num,nZuc = readlattice(lfile_z,numstart=0)
	ZG,ZG_num,nZGuc = readlattice(lfile_ig,numstart=nZuc)
	Gas,Gas_num,nGuc = readlattice(lfile_g,numstart=nZuc+nZGuc)
	print nZuc,nZGuc,nGuc
	label = []
	if (nZuc > 0):
		row = 0
		while row < nZuc:
			label+=["ZI"]*15
			label+=["ZZ","ZZ"]
			label+=["ZS","ZS","ZS","ZS","ZS"]
			label+=["ZZ","ZZ"]
			label+=["ZS","ZS","ZS","ZS","ZS"]
			label+=["ZZ","ZZ"]
			row+=+1
	if (nZGuc > 0):
		row = 0
		while row < nZGuc:
			label+=["ZGI"]*15
			row+=1
	if (nGuc > 0):
		row = 0
		while row < nGuc:
			label+=["GI"]*15
			row+=1
	print np.unique(label)
	np.save(sys.argv[1]+'_label',label)
	del label
	
	rows = 31*nZuc+15*nZGuc+15*nGuc
	neighbors = np.memmap('testmmp_2.npy',dtype=int,mode='w+',shape=(rows,10))
	print "empty neighbors created"
	print "neighbors initialized with -1"
	i = 1
	while i <= rows:
		neighbors[i-1,0] = i
		i = i+1
	num_neighbors = 0*neighbors[:,0]
	print "neighbors initilized with index"

	if (nZuc > 0):
		connection = {}
		connection[0]=[15,17]
		connection[1]=[18]
		connection[2]=[19]
		connection[3]=[16,20]
		connection[4]=[15,16,21]
		connection[5]=[17,24]
		connection[6]=[18,22,25]
		connection[7]=[19,23,26]
		connection[8]=[20,27]
		connection[9]=[21,22,23,28]
		connection[10]=[24,29]
		connection[11]=[25]
		connection[12]=[26]
		connection[13]=[27,30]
		connection[14]=[28,29,30]
		#########################
		connection[15]=[0,4]
		connection[16]=[3,4]
		connection[17]=[0,5]
		connection[18]=[1,6]
		connection[19]=[2,7]
		connection[20]=[3,8]
		connection[21]=[4,9]
		connection[22]=[6,9]
		connection[23]=[7,9]
		connection[24]=[5,10]
		connection[25]=[6,11]
		connection[26]=[7,12]
		connection[27]=[8,13]
		connection[28]=[9,14]
		connection[29]=[10,14]
		connection[30]=[13,14]
		for key in connection:
			for value in connection[key]:
				num_neighbors[int(key):31*nZuc:31]+=1
				neighbors[int(key):31*nZuc:31,num_neighbors[int(key)]]=neighbors[int(value):31*nZuc:31,0]
	print "Zeolite part initialized"
	if (nZGuc > 0):
		connection = {}
		connection[0]=[4,5]
		connection[1]=[4,6]
		connection[2]=[4,7]
		connection[3]=[4,8]
		connection[4]=[0,1,2,3,9]
		connection[5]=[0,10,9]
		connection[6]=[1,11,9]
		connection[7]=[2,12,9]
		connection[8]=[3,9,13]
		connection[9]=[4,5,6,7,8,14]
		connection[10]=[5,14]
		connection[11]=[14,6]
		connection[12]=[7,14]
		connection[13]=[8,14]
		connection[14]=[9,10,11,12,13]
		for key in connection:
			for value in connection[key]:
				k = int(key)+31*nZuc
				v = int(value)+31*nZuc
				upto = 31*nZuc+15*nZGuc
				num_neighbors[k:upto:15]+=1
				neighbors[k:upto:15,num_neighbors[k]]=neighbors[v:upto:15,0]
	print "Interface part initialized"
	if (nGuc > 0):
		connection = {}
		connection[0]=[4,5]
		connection[1]=[4,6]
		connection[2]=[4,7]
		connection[3]=[4,8]
		connection[4]=[0,1,2,3,9]
		connection[5]=[0,10,9]
		connection[6]=[1,11,9]
		connection[7]=[2,12,9]
		connection[8]=[3,9,13]
		connection[9]=[4,5,6,7,8,14]
		connection[10]=[5,14]
		connection[11]=[14,6]
		connection[12]=[7,14]
		connection[13]=[8,14]
		connection[14]=[9,10,11,12,13]
		for key in connection:
			for value in connection[key]:
				k = int(key)+31*nZuc+15*nZGuc
				v = int(value)+31*nZuc+15*nZGuc
				upto = 31*nZuc+15*nZGuc+15*nGuc
				num_neighbors[k:upto:15]+=1
				neighbors[k:upto:15,num_neighbors[k]]=neighbors[v:upto:15,0]
	print "Gas part initialized"
	print neighbors
	print np.unique(num_neighbors)
	del connection


# In[7]:


def getid(num):
	global nZuc,nZGuc,nGuc
	if (num <= nZuc):
		id = 31*(num-1)
	elif (num <= nZuc+nZGuc):
		id = 31*nZuc + 15*(num-nZuc-1)
	elif (num <= nZuc+nZGuc+nGuc):
		id = 31*nZuc + 15*(nZGuc) + 15*(num-nZuc-nZGuc-1)
	return id


# In[8]:


def modify(l1,l2,id0,id1):
	global neighbors
	for a,b in zip(l1,l2):
		n1 = neighbors[id0+a-1,0]
		n2 = neighbors[id1+b-1,0]
		neighbors[id0+a-1,0] = min(n1,n2)
		neighbors[id1+b-1,0] = min(n1,n2)


# In[9]:


def periodic_neighborlist_g():
	global siteID, neighbors, Gas, Gas_num
	arr = np.load(sys.argv[1]+'_boundary.npy')
	xyz = Zeo+ZG+Gas
	xyz_num = Zeo_num+ZG_num+Gas_num
	Lx = np.size(Gas,axis=0)
	Ly = np.size(Gas,axis=1)
	Lz = np.size(Gas,axis=2)
	print Lx,Ly,Lz
	for i,j,k in zip(*np.where(arr==1)):
		if (Gas[i,j,k] == 1):
			n = int(Gas_num[i,j,k])
			id0 = getid(Gas_num[i,j,k])
			#mask = siteID[n-1]-1
			i1 = i-1
			i2 = i+1
			j1 = j-1
			j2 = j+1
			k1 = k-1
			k2 = k+1
			if (k1<0):
				px = i
				py = j
				pz = Lz-1
				a = [6,9,11,14,1,4]
				b = [10,10,15,15,5,5]
				id1 = getid(Gas_num[px,py,pz])
				for m,l in zip(a,b):
					loc=np.where(neighbors[id0+m-1:id0+m,1:]==neighbors[id1+l-1,0])[1]
					loc2 = np.where(neighbors[id0+m-1:id0+m,1:]==0)[1]
					if (~len(loc)):
						neighbors[id0+m-1,loc2[0]+1]=neighbors[id1+l-1,0]
						num_neighbors[id0+m-1]+=1
				#if (Gas[px,py,pz]==1):
				#	for l,m in zip(a,b):
				#		if (siteID[Gas_num[px,py,pz]-1][m-1] not in neighbors[mask[l-1]][1]):
				#			neighbors[mask[l-1]][1]+=[siteID[Gas_num[px,py,pz]-1][m-1]]
			if (i1<0):
				px = Lx-1
				py = j
				pz = k
				a = [11,7,1,12,6,2]
				b = [15,10,5,15,10,5]
				id1 = getid(Gas_num[px,py,pz])
				for m,l in zip(a,b):
					loc=np.where(neighbors[id0+m-1:id0+m,1:]==neighbors[id1+l-1,0])[1]
					loc2 = np.where(neighbors[id0+m-1:id0+m,1:]==0)[1]
					if (~len(loc)):
						neighbors[id0+m-1,loc2[0]+1]=neighbors[id1+l-1,0]
						num_neighbors[id0+m-1]+=1

				#if (Gas[px,py,pz]==1):
				#	for l,m in zip(a,b):
				#		if (siteID[Gas_num[px,py,pz]-1][m-1] not in neighbors[mask[l-1]][1]):
				#			neighbors[mask[l-1]][1]+=[siteID[Gas_num[px,py,pz]-1][m-1]]
			if (j1<0):
				px = i
				py = Ly-1
				pz = k
				a = [1,2,5,3,4]
				b = [6,7,10,8,9]
				id1 = getid(Gas_num[px,py,pz])
				for m,l in zip(a,b):
					loc=np.where(neighbors[id0+m-1:id0+m,1:]==neighbors[id1+l-1,0])[1]
					loc2 = np.where(neighbors[id0+m-1:id0+m,1:]==0)[1]
					if (~len(loc)):
						neighbors[id0+m-1,loc2[0]+1]=neighbors[id1+l-1,0]
						num_neighbors[id0+m-1]+=1

				#if (Gas[px,py,pz]==1):
				#	for l,m in zip(a,b):
				#		if (siteID[Gas_num[px,py,pz]-1][m-1] not in neighbors[mask[l-1]][1]):
				#			neighbors[mask[l-1]][1]+=[siteID[Gas_num[px,py,pz]-1][m-1]]
			if (i2>=Lx):
				px = 0
				py = j
				pz = k
				a = [14,8,4,13,9,3]
				b = [15,10,5,15,10,5]
				id1 = getid(Gas_num[px,py,pz])
				for m,l in zip(a,b):
					loc=np.where(neighbors[id0+m-1:id0+m,1:]==neighbors[id1+l-1,0])[1]
					loc2 = np.where(neighbors[id0+m-1:id0+m,1:]==0)[1]
					if (~len(loc)):
						neighbors[id0+m-1,loc2[0]+1]=neighbors[id1+l-1,0]
						num_neighbors[id0+m-1]+=1

				#if (Gas[px,py,pz]==1):
				#	for l,m in zip(a,b):
				#		if (siteID[Gas_num[px,py,pz]-1][m-1] not in neighbors[mask[l-1]][1]):
				#			neighbors[mask[l-1]][1]+=[siteID[Gas_num[px,py,pz]-1][m-1]]
			if (j2>=Ly):
				px = i
				py = 0
				pz = k
				a = [11,12,15,13,14]
				b = [6,7,10,8,9]
				id1 = getid(Gas_num[px,py,pz])
				for m,l in zip(a,b):
					loc=np.where(neighbors[id0+m-1:id0+m,1:]==neighbors[id1+l-1,0])[1]
					loc2 = np.where(neighbors[id0+m-1:id0+m,1:]==0)[1]
					if (~len(loc)):
						neighbors[id0+m-1,loc2[0]+1]=neighbors[id1+l-1,0]
						num_neighbors[id0+m-1]+=1

				#if (Gas[px,py,pz]==1):
				#	for l,m in zip(a,b):
				#		if (siteID[Gas_num[px,py,pz]-1][m-1] not in neighbors[mask[l-1]][1]):
				#			neighbors[mask[l-1]][1]+=[siteID[Gas_num[px,py,pz]-1][m-1]]
			if (k2>=Lz):
				px = i
				py = j
				pz = 0
				a = [12,2,13,3,7,8]
				b = [15,5,15,5,10,10]
				id1 = getid(Gas_num[px,py,pz])
				for m,l in zip(a,b):
					loc=np.where(neighbors[id0+m-1:id0+m,1:]==neighbors[id1+l-1,0])[1]
					loc2 = np.where(neighbors[id0+m-1:id0+m,1:]==0)[1]
					if (~len(loc)):
						neighbors[id0+m-1,loc2[0]+1]=neighbors[id1+l-1,0]
						num_neighbors[id0+m-1]+=1

				#if (Gas[px,py,pz]==1):
				#	for l,m in zip(a,b):
				#		if (siteID[Gas_num[px,py,pz]-1][m-1] not in neighbors[mask[l-1]][1]):
				#			neighbors[mask[l-1]][1]+=[siteID[Gas_num[px,py,pz]-1][m-1]]
			if (i1<0 and k1<0):
				px = Lx-1
				py = j
				pz = Lz-1
				a = [6,11,1]
				b = [10,15,5]
				id1 = getid(Gas_num[px,py,pz])
				for m,l in zip(a,b):
					loc=np.where(neighbors[id0+m-1:id0+m,1:]==neighbors[id1+l-1,0])[1]
					loc2 = np.where(neighbors[id0+m-1:id0+m,1:]==0)[1]
					if (~len(loc)):
						neighbors[id0+m-1,loc2[0]+1]=neighbors[id1+l-1,0]
						num_neighbors[id0+m-1]+=1

					#if (Gas[px,py,pz]==1):
					#	for l,m in zip(a,b):
					#		if (siteID[Gas_num[px,py,pz]-1][m-1] not in neighbors[mask[l-1]][1]):
					#			neighbors[mask[l-1]][1]+=[siteID[Gas_num[px,py,pz]-1][m-1]]
			if (i2>=Lx and k2>=Lz):
				px = 0
				py = j
				pz = 0
				a = [13,3,8]
				b = [15,5,10]
				id1 = getid(Gas_num[px,py,pz])
				for m,l in zip(a,b):
					loc=np.where(neighbors[id0+m-1:id0+m,1:]==neighbors[id1+l-1,0])[1]
					loc2 = np.where(neighbors[id0+m-1:id0+m,1:]==0)[1]
					if (~len(loc)):
						neighbors[id0+m-1,loc2[0]+1]=neighbors[id1+l-1,0]
						num_neighbors[id0+m-1]+=1

					#if (Gas[px,py,pz]==1):
					#	for l,m in zip(a,b):
					#		if (siteID[Gas_num[px,py,pz]-1][m-1] not in neighbors[mask[l-1]][1]):
					#			neighbors[mask[l-1]][1]+=[siteID[Gas_num[px,py,pz]-1][m-1]]
			if (i2>=Lx and k1<0):
				px = 0
				py = j
				pz = Lz-1
				a = [9,14,4]
				b = [10,15,5]
				id1 = getid(Gas_num[px,py,pz])
				for m,l in zip(a,b):
					loc=np.where(neighbors[id0+m-1:id0+m,1:]==neighbors[id1+l-1,0])[1]
					loc2 = np.where(neighbors[id0+m-1:id0+m,1:]==0)[1]
					if (~len(loc)):
						neighbors[id0+m-1,loc2[0]+1]=neighbors[id1+l-1,0]
						num_neighbors[id0+m-1]+=1

					#if (Gas[px,py,pz]==1):
					#	for l,m in zip(a,b):
					#		if (siteID[Gas_num[px,py,pz]-1][m-1] not in neighbors[mask[l-1]][1]):
					#			neighbors[mask[l-1]][1]+=[siteID[Gas_num[px,py,pz]-1][m-1]]
			if (i1<0 and k2>=Lz):
				px = Lx-1
				py = j
				pz = 0
				a = [12,2,7]
				b = [15,5,10]
				id1 = getid(Gas_num[px,py,pz])
				for m,l in zip(a,b):
					loc=np.where(neighbors[id0+m-1:id0+m,1:]==neighbors[id1+l-1,0])[1]
					loc2 = np.where(neighbors[id0+m-1:id0+m,1:]==0)[1]
					if (~len(loc)):
						neighbors[id0+m-1,loc2[0]+1]=neighbors[id1+l-1,0]
						num_neighbors[id0+m-1]+=1
					#if (Gas[px,py,pz]==1):
					#	for l,m in zip(a,b):
					#		if (siteID[Gas_num[px,py,pz]-1][m-1] not in neighbors[mask[l-1]][1]):
					#			neighbors[mask[l-1]][1]+=[siteID[Gas_num[px,py,pz]-1][m-1]]
			#del mask
			#gc.collect()
	del arr
	#return alist,nlist


# In[10]:


def periodic_neighborlist_z(alist,nlist,xyz,xyz_num):
	Lx = np.size(xyz,axis=0)
	Ly = np.size(xyz,axis=1)
	Lz = np.size(xyz,axis=2)
	print Lx,Ly,Lz
	for i in range(0,Lx,1):
		for j in range(0,Ly,1):
			for k in range(0,Lz,1):
				if (xyz[i,j,k] == 1):
					n = int(xyz_num[i,j,k])
					mask = alist[n-1]-1
					i1 = i-1
					i2 = i+1
					j1 = j-1
					j2 = j+1
					k1 = k-1
					k2 = k+1
					if (k1<0):
						px = i
						py = j
						pz = Lz-1
						a = [6,9]
						b = [23,24]
						if (xyz[px,py,pz]==1):
							for l,m in zip(a,b):
								if (alist[xyz_num[px,py,pz]-1][m-1] not in nlist[mask[l-1]][1]):
									nlist[mask[l-1]][1]+=[alist[xyz_num[px,py,pz]-1][m-1]]
					if (i1<0):
						px = Lx-1
						py = j
						pz = k
						a = [11,7,1]
						b = [31,24,17]
						if (xyz[px,py,pz]==1):
							for l,m in zip(a,b):
								if (alist[xyz_num[px,py,pz]-1][m-1] not in nlist[mask[l-1]][1]):
									nlist[mask[l-1]][1]+=[alist[xyz_num[px,py,pz]-1][m-1]]
					if (j1<0):
						px = i
						py = Ly-1
						pz = k
						a = [1,2,5,3,4]
						b = [25,26,29,27,28]
						if (xyz[px,py,pz]==1):
							for l,m in zip(a,b):
								if (alist[xyz_num[px,py,pz]-1][m-1] not in nlist[mask[l-1]][1]):
									nlist[mask[l-1]][1]+=[alist[xyz_num[px,py,pz]-1][m-1]]
					if (i2>=Lx):
						px = 0
						py = j
						pz = k
						a = [14,8,4]
						b = [30,23,16]
						if (xyz[px,py,pz]==1):
							for l,m in zip(a,b):
								if (alist[xyz_num[px,py,pz]-1][m-1] not in nlist[mask[l-1]][1]):
									nlist[mask[l-1]][1]+=[alist[xyz_num[px,py,pz]-1][m-1]]
					if (j2>=Ly):
						px = i
						py = 0
						pz = k
						a = [11,12,15,13,14]
						b = [18,19,22,20,21]
						if (xyz[px,py,pz]==1):
							for l,m in zip(a,b):
								if (alist[xyz_num[px,py,pz]-1][m-1] not in nlist[mask[l-1]][1]):
									nlist[mask[l-1]][1]+=[alist[xyz_num[px,py,pz]-1][m-1]]
					if (k2>=Lz):
						px = i
						py = j
						pz = 0
						a = [12,2,13,3]
						b = [30,16,31,17]
						if (xyz[px,py,pz]==1):
							for l,m in zip(a,b):
								if (alist[xyz_num[px,py,pz]-1][m-1] not in nlist[mask[l-1]][1]):
									nlist[mask[l-1]][1]+=[alist[xyz_num[px,py,pz]-1][m-1]]
					if (i1<0 and k1<0):
						px = Lx-1
						py = j
						pz = Lz-1
						a = [6]
						b = [24]
						if (xyz[px,py,pz]==1):
							for l,m in zip(a,b):
								if (alist[xyz_num[px,py,pz]-1][m-1] not in nlist[mask[l-1]][1]):
									nlist[mask[l-1]][1]+=[alist[xyz_num[px,py,pz]-1][m-1]]
					if (i2>=Lx and k2>=Lz):
						px = 0
						py = j
						pz = 0
						a = [13,3]
						b = [30,16]
						if (xyz[px,py,pz]==1):
							for l,m in zip(a,b):
								if (alist[xyz_num[px,py,pz]-1][m-1] not in nlist[mask[l-1]][1]):
									nlist[mask[l-1]][1]+=[alist[xyz_num[px,py,pz]-1][m-1]]
					if (i2>=Lx and k1<0):
						px = 0
						py = j
						pz = Lz-1
						a = [9]
						b = [23]
						if (xyz[px,py,pz]==1):
							for l,m in zip(a,b):
								if (alist[xyz_num[px,py,pz]-1][m-1] not in nlist[mask[l-1]][1]):
									nlist[mask[l-1]][1]+=[alist[xyz_num[px,py,pz]-1][m-1]]
					if (i1<0 and k2>=Lz):
						px = Lx-1
						py = j
						pz = 0
						a = [12,2]
						b = [31,17]
						if (xyz[px,py,pz]==1):
							for l,m in zip(a,b):
								if (alist[xyz_num[px,py,pz]-1][m-1] not in nlist[mask[l-1]][1]):
									nlist[mask[l-1]][1]+=[alist[xyz_num[px,py,pz]-1][m-1]]
	return alist,nlist


# In[13]:


def remove_neighborlist_if_z():
	global siteID,neighbors,Zeo,Zeo_num,ZG,ZG_num,nZuc,nZGuc,nGuc,num_neighbors
	Lx = np.size(Zeo,axis=0)
	Ly = np.size(Zeo,axis=1)
	Lz = np.size(Zeo,axis=2)
	for i in range(0,Lx,1):
		for j in range(0,Ly,1):
			for k in range(0,Lz,1):
				if (ZG[i,j,k] == 1):
					n = int(ZG_num[i,j,k])
					id = int(31*nZuc+15*(ZG_num[i,j,k]-nZuc-1))
					i1 = i-1
					i2 = i+1
					j1 = j-1
					j2 = j+1
					k1 = k-1
					k2 = k+1
					if (i1>=0):
						if (Zeo[i1,j,k] == 1):
							a = [11,6,1,12,7,2,15,15,10,10,5,5,11,6,1,6,12,7,7,2]
							b = [15,10,5,15,10,5,11,12,6,7,1,2,6,11,6,1,7,12,2,7]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1
					if (i2<np.size(Zeo,axis=0)):
						if (Zeo[i2,j,k] == 1):
							a = [15,10,5,15,10,5,14,9,4,13,8,3,14,9,9,4,13,8,8,3]
							b = [14,9,4,13,8,3,15,10,5,15,10,5,9,14,4,9,8,13,3,8]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1
					if (j1>=0):
						if (Zeo[i,j1,k] == 1):
							a = [5,2,5,3,5,1,5,4]
							b = [2,5,3,5,1,5,4,5]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1
					if (j2<np.size(Zeo,axis=1)):
						if (Zeo[i,j2,k] == 1):
							a = [15,12,15,13,15,11,15,14]
							b = [12,15,13,15,11,15,14,15]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1
					if (k1>=0):
						if (Zeo[i,j,k1] == 1):
							a = [11,15,14,15,6,10,9,10,1,5,4,5,11,6,6,1,14,9,9,4]
							b = [15,11,15,14,10,6,10,9,5,1,5,4,6,11,1,6,9,14,4,9]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1
					if (k2<np.size(Zeo,axis=2)):
						if (Zeo[i,j,k2] == 1):
							a = [12,15,13,15,7,10,8,10,2,5,3,5,12,7,7,2,13,8,8,3]
							b = [15,12,15,13,10,7,10,8,5,2,5,3,7,12,2,7,8,13,3,8]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1
					if (i1>=0 and j1>=0):
						if (Zeo[i1,j1,k] == 1):
							a = [5,2,1,5]
							b = [2,5,5,1]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1
					if (i2<np.size(Zeo,axis=0) and j1>=0):
						if (Zeo[i2,j1,k] == 1):
							a = [5,4,5,3]
							b = [4,5,3,5]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (i1>=0 and j2<np.size(Zeo,axis=1)):
						if (Zeo[i1,j2,k] == 1):
							a = [11,15,12,15]
							b = [15,11,15,12]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (i2<np.size(Zeo,axis=0) and j2<np.size(Zeo,axis=1)):
						if (Zeo[i2,j2,k] == 1):
							a = [15,13,15,14]
							b = [13,15,14,15]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (i1>=0 and k1>=0):
						if (Zeo[i1,j,k1] == 1):
							a = [11,15,1,5,11,6,6,1,6,10]
							b = [15,11,5,1,6,11,1,6,10,6]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (i2<np.size(Zeo,axis=0) and k1>=0):
						if (Zeo[i2,j,k1] == 1):
							a = [15,14,5,4,14,9,9,4,10,9]
							b = [14,15,4,5,9,14,4,9,9,10]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (i1>=0 and k2<np.size(Zeo,axis=2)):
						if (Zeo[i1,j,k2] == 1):
							a = [7,10,12,7,7,2,12,15,2,5]
							b = [10,7,7,12,2,7,15,12,5,2]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (i2<np.size(Zeo,axis=0) and k2<np.size(Zeo,axis=2)):
						if (Zeo[i2,j,k2] == 1):
							a = [10,8,13,8,8,3,15,13,5,3]
							b = [8,10,8,13,3,8,13,15,3,5]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (j1>=0 and k1>=0):
						if (Zeo[i,j1,k1] == 1):
							a = [1,5,4,5]
							b = [5,1,5,4]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (j2<np.size(Zeo,axis=1) and k1>=0):
						if (Zeo[i,j2,k1] == 1):
							a = [11,15,14,15]
							b = [15,11,15,14]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1
					if (j1>=0 and k2<np.size(Zeo,axis=2)):
						if (Zeo[i,j1,k2] == 1):
							a = [2,5,3,5]
							b = [5,2,5,3]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (j2<np.size(Zeo,axis=1) and k2<np.size(Zeo,axis=2)):
						if (Zeo[i,j2,k2] == 1):
							a = [12,15,13,15]
							b = [15,12,15,13]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (i1>=0 and j1>=0 and k1>=0):
						if (Zeo[i1,j1,k1] == 1):
							a = [1,5]
							b = [5,1]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (i2<np.size(Zeo,axis=0) and j1>=0 and k1>=0):
						if (Zeo[i2,j1,k1] == 1):
							a = [4,5]
							b = [5,4]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (i1>=0 and j2<np.size(Zeo,axis=1) and k1>=0):
						if (Zeo[i1,j2,k1] == 1):
							a = [11,15]
							b = [15,11]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (i1>=0 and j1>=0 and k2<np.size(Zeo,axis=2)):
						if (Zeo[i1,j1,k2] == 1):
							a = [2,5]
							b = [5,2]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (i2<np.size(Zeo,axis=0) and j2<np.size(Zeo,axis=1) and k1>=0):
						if (Zeo[i2,j2,k1] == 1):
							a = [14,15]
							b = [15,14]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (i1>=0 and j2<np.size(Zeo,axis=1) and k2<np.size(Zeo,axis=2)):
						if (Zeo[i1,j2,k2] == 1):
							a = [12,15]
							b = [15,12]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (i2<np.size(Zeo,axis=0) and j1>=0 and k2<np.size(Zeo,axis=2)):
						if (Zeo[i2,j1,k2] == 1):
							a = [5,3]
							b = [3,5]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1

					if (i2<np.size(Zeo,axis=0) and j2<np.size(Zeo,axis=1) and k2<np.size(Zeo,axis=2)):
						if (Zeo[i2,j2,k2] == 1):
							a = [15,13]
							b = [13,15]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								if (len(loc)):
									neighbors[id+m-1,loc[0]+1]=0
									num_neighbors[id+m-1]-=1
	#del i,j,k,i1,i2,j1,j2,k1,k2,a,b,Lx,Ly,Lz
	#return alist,nlist


# In[14]:


def add_neighborlist_if_z():
	global siteID,neighbors,Zeo,Zeo_num,ZG,ZG_num,nZuc,nZGuc,nGuc,num_neighbors
	Lx = np.size(Zeo,axis=0)
	Ly = np.size(Zeo,axis=1)
	Lz = np.size(Zeo,axis=2)
	for i in range(0,Lx,1):
		for j in range(0,Ly,1):
			for k in range(0,Lz,1):
				if (ZG[i,j,k] == 1):
					n = int(ZG_num[i,j,k])
					id = int(31*nZuc+15*(ZG_num[i,j,k]-nZuc-1))
					i1 = i-1
					i2 = i+1
					j1 = j-1
					j2 = j+1
					k1 = k-1
					k2 = k+1
					if (i1>=0 and k1>=0):
						if (Zeo[i1,j,k1] == 1):
							a = [6,10]
							b = [10,6]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								loc2 = np.where(neighbors[id+m-1:id+m,1:]==0)[1]
								if (~len(loc)):
									neighbors[id+m-1,loc2[0]+1]=neighbors[id+l-1,0]
									num_neighbors[id+m-1]+=1

					if (i2<np.size(Zeo,axis=0) and k1>=0):
						if (Zeo[i2,j,k1] == 1):
							a = [10,9]
							b = [9,10]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								loc2 = np.where(neighbors[id+m-1:id+m,1:]==0)[1]
								if (~len(loc)):
									neighbors[id+m-1,loc2[0]+1]=neighbors[id+l-1,0]
									num_neighbors[id+m-1]+=1
	
					if (i1>=0 and k2<np.size(Zeo,axis=2)):
						if (Zeo[i1,j,k2] == 1):
							a = [12,15,2,5]
							b = [15,12,5,2]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								loc2 = np.where(neighbors[id+m-1:id+m,1:]==0)[1]
								if (~len(loc)):
									neighbors[id+m-1,loc2[0]+1]=neighbors[id+l-1,0]
									num_neighbors[id+m-1]+=1

					if (i2<np.size(Zeo,axis=0) and k2<np.size(Zeo,axis=2)):
						if (Zeo[i2,j,k2] == 1):
							a = [15,13,5,3]
							b = [13,15,3,5]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								loc2 = np.where(neighbors[id+m-1:id+m,1:]==0)[1]
								if (~len(loc)):
									neighbors[id+m-1,loc2[0]+1]=neighbors[id+l-1,0]
									num_neighbors[id+m-1]+=1

					if (i1>=0 and j1>=0 and k2<np.size(Zeo,axis=2)):
						if (Zeo[i1,j1,k2] == 1):
							a = [2,5]
							b = [5,2]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								loc2 = np.where(neighbors[id+m-1:id+m,1:]==0)[1]
								if (~len(loc)):
									neighbors[id+m-1,loc2[0]+1]=neighbors[id+l-1,0]
									num_neighbors[id+m-1]+=1

					if (i1>=0 and j2<np.size(Zeo,axis=1) and k2<np.size(Zeo,axis=2)):
						if (Zeo[i1,j2,k2] == 1):
							a = [12,15]
							b = [15,12]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								loc2 = np.where(neighbors[id+m-1:id+m,1:]==0)[1]
								if (~len(loc)):
									neighbors[id+m-1,loc2[0]+1]=neighbors[id+l-1,0]
									num_neighbors[id+m-1]+=1

					if (i2<np.size(Zeo,axis=0) and j1>=0 and k2<np.size(Zeo,axis=2)):
						if (Zeo[i2,j1,k2] == 1):
							a = [3,5]
							b = [5,3]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								loc2 = np.where(neighbors[id+m-1:id+m,1:]==0)[1]
								if (~len(loc)):
									neighbors[id+m-1,loc2[0]+1]=neighbors[id+l-1,0]
									num_neighbors[id+m-1]+=1

					if (i2<np.size(Zeo,axis=0) and j2<np.size(Zeo,axis=1) and k2<np.size(Zeo,axis=2)):
						if (Zeo[i2,j2,k2] == 1):
							a = [13,15]
							b = [15,13]
							for l,m in zip(a,b):
								loc=np.where(neighbors[id+m-1:id+m,1:]==neighbors[id+l-1,0])[1]
								loc2 = np.where(neighbors[id+m-1:id+m,1:]==0)[1]
								if (~len(loc)):
									neighbors[id+m-1,loc2[0]+1]=neighbors[id+l-1,0]
									num_neighbors[id+m-1]+=1
	#del i,j,k,i1,i2,j1,j2,k1,k2,a,b,Lx,Ly,Lz
	#return alist,nlist


# In[15]:


def nonperiodic_neighborlist_2(): #xyz,xyz_num,Z):
	global siteID,neighbors,Zeo,ZG,Gas,ZG_num,Zeo_num,Gas_num,nZuc,nZGuc,nGuc
	xyz = Zeo+ZG+Gas
	xyz_num = Zeo_num+ZG_num+Gas_num
	print xyz
	print xyz_num
	print nZuc,nZGuc,nGuc
	Lx = np.size(xyz,axis=0)
	Ly = np.size(xyz,axis=1)
	Lz = np.size(xyz,axis=2)
	z = 0
	print "Started nonperiodic masking"
	print "Created empty mas array"
	i = 0
	while i<Lx:
	#for i in range(0,Lx,1):
		#print i
		#start = time.time()
		j = 0
		while j<Ly:
		#for j in range(0,Ly,1):
			k = 0
			while k<Lz:
			#for k in range(0,Lz,1):
				#n = int(xyz_num[i,j,k])
				id0 = getid(xyz_num[i,j,k])
				#print xyz_num[i,j,k]
				#mask = neighbors[id0+1-1:idf,0]
				i1 = i-1
				i2 = i+1
				j1 = j-1
				j2 = j+1
				k1 = k-1
				k2 = k+1
				if (i1>=0):
					l1 = [12,11,6,7,1,2]
					l2 = [13,14,9,8,4,3]
					id1 = getid(xyz_num[i1,j,k])
					modify(l1,l2,id0,id1)
					if(Zeo[i,j,k]*Zeo[i1,j,k]==1):
						l1 = [25,26,18,19]
						l2 = [28,27,21,20]
						modify(l1,l2,id0,id1)
				if (i2<Lx):
					l1 = [13,14,9,8,4,3]
					l2 = [12,11,6,7,1,2]
					id1 = getid(xyz_num[i2,j,k])
					modify(l1,l2,id0,id1)
					if (Zeo[i,j,k]*Zeo[i2,j,k] == 1):
						l2 = [25,26,18,19]
						l1 = [28,27,21,20]
						modify(l1,l2,id0,id1)
				if (j1>=0):
					l1 = [1,2,5,4,3]
					l2 = [11,12,15,14,13]
					id1 = getid(xyz_num[i,j1,k])
					modify(l1,l2,id0,id1)
					if (Zeo[i,j,k]*Zeo[i,j1,k] == 1):
						l1 = [16,17]
						l2 = [30,31]
						modify(l1,l2,id0,id1)
				if (j2<Ly):
					l1 = [11,12,15,14,13]
					l2 = [1,2,5,4,3]
					id1 = getid(xyz_num[i,j2,k])
					modify(l1,l2,id0,id1)
					if (Zeo[i,j,k]*Zeo[i,j2,k] == 1):
						l1 = [30,31]
						l2 = [16,17]
						modify(l1,l2,id0,id1)
				if (k1>=0):
					l1 = [11,14,6,9,1,4]
					l2 = [12,13,7,8,2,3]
					id1 = getid(xyz_num[i,j,k1])
					modify(l1,l2,id0,id1)
					if (Zeo[i,j,k]*Zeo[i,j,k1] == 1):
						l1 = [25,28,18,21]
						l2 = [26,27,19,20]
						modify(l1,l2,id0,id1)
				if (k2<Lz):
					l1 = [12,13,7,8,2,3]
					l2 = [11,14,6,9,1,4]
					id1 = getid(xyz_num[i,j,k2])
					modify(l1,l2,id0,id1)
					if (Zeo[i,j,k]*Zeo[i,j,k2] == 1):
						l1 = [26,27,19,20]
						l2 = [25,28,18,21]
						modify(l1,l2,id0,id1)
				if (i1>=0 and j1>=0):
					l1 = [1,2]
					l2 = [14,13]
					id1 = getid(xyz_num[i1,j1,k])
					modify(l1,l2,id0,id1)
				if (i2<Lx and j1>=0):
					l1 = [3,4]
					l2 = [12,11]
					id1 = getid(xyz_num[i2,j1,k])
					modify(l1,l2,id0,id1)
				if (i1>=0 and j2<Ly):
					l1 = [11,12]
					l2 = [4,3]
					id1 = getid(xyz_num[i1,j2,k])
					modify(l1,l2,id0,id1)
				if (i2<Lx and j2<Ly):
					l1 = [14,13]
					l2 = [1,2]
					id1 = getid(xyz_num[i2,j2,k])
					modify(l1,l2,id0,id1)
				if (i1>=0 and k1>=0):
					l1 = [11,6,1]
					l2 = [13,8,3]
					id1 = getid(xyz_num[i1,j,k1])
					modify(l1,l2,id0,id1)
					if (Zeo[i,j,k]*Zeo[i1,j,k1] == 1):
						l1 = [25,18]
						l2 = [27,20]
						modify(l1,l2,id0,id1)
				if (i2<Lx and k1>=0):
					l1 = [14,9,4]
					l2 = [12,7,2]
					id1 = getid(xyz_num[i2,j,k1])
					modify(l1,l2,id0,id1)
					if (Zeo[i,j,k]*Zeo[i2,j,k1] == 1):
						l1 = [28,21]
						l2 = [26,19]
						modify(l1,l2,id0,id1)
				if (i1>=0 and k2<Lz):
					l1 = [12,7,2]
					l2 = [14,9,4]
					id1 = getid(xyz_num[i1,j,k2])
					modify(l1,l2,id0,id1)
					if (Zeo[i,j,k]*Zeo[i1,j,k2] == 1):
						l1 = [26,19]
						l2 = [28,21]
						modify(l1,l2,id0,id1)
				if (i2<Lx and k2<Lz):
					l1 = [13,8,3]
					l2 = [11,6,1]
					id1 = getid(xyz_num[i2,j,k2])
					modify(l1,l2,id0,id1)
					if (Zeo[i,j,k]*Zeo[i2,j,k2] == 1):
						l1 = [27,20]
						l2 = [25,18]
						modify(l1,l2,id0,id1)
				if (j1>=0 and k1>=0):
					l1 = [1,4]
					l2 = [12,13]
					id1 = getid(xyz_num[i,j1,k1])
					modify(l1,l2,id0,id1)
				if (j2<Ly and k1>=0):
					l1 = [11,14]
					l2 = [2,3]
					id1 = getid(xyz_num[i,j2,k1])
					modify(l1,l2,id0,id1)
				if (j1>=0 and k2<Lz):
					l1 = [2,3]
					l2 = [11,14]
					id1 = getid(xyz_num[i,j1,k2])
					modify(l1,l2,id0,id1)
				if (j2<Ly and k2<Lz):
					l1 = [12,13]
					l2 = [1,4]
					id1 = getid(xyz_num[i,j2,k2])
					modify(l1,l2,id0,id1)
				if (i1>=0 and j1>=0 and k1>=0):
					id1 = getid(xyz_num[i1,j1,k1])
					modify([1],[13],id0,id1)
				if (i2<Lx and j1>=0 and k1>=0):
					id1 = getid(xyz_num[i2,j1,k1])
					modify([4],[12],id0,id1)
				if (i1>=0 and j2<Ly and k1>=0):
					id1 = getid(xyz_num[i1,j2,k1])
					modify([11],[3],id0,id1)
				if (i1>=0 and j1>=0 and k2<Lz):
					id1 = getid(xyz_num[i1,j1,k2])
					modify([2],[14],id0,id1)
				if (i2<Lx and j2<Ly and k1>=0):
					id1 = getid(xyz_num[i2,j2,k1])
					modify([14],[2],id0,id1)
				if (i1>=0 and j2<Ly and k2<Lz):
					id1 = getid(xyz_num[i1,j2,k2])
					modify([12],[4],id0,id1)
				if (i2<Lx and j1>=0 and k2<Lz):
					id1 = getid(xyz_num[i2,j1,k2])
					modify([3],[11],id0,id1)
				if (i2<Lx and j2<Ly and k2<Lz):
					id1 = getid(xyz_num[i2,j2,k2])
					modify([13],[1],id0,id1)
				k = k+1
			j = j+1
		i = i+1
		#print Lx-i,time.time()#-start
	#return alist,nlist


# In[16]:


def nonperiodic_neighborlist(alist,nlist,xyz,xyz_num,alist2):
	Lx = np.size(xyz,axis=0)
	Ly = np.size(xyz,axis=1)
	Lz = np.size(xyz,axis=2)
	for i,j,k in zip(*np.where(xyz==1)):#range(0,Lx,1):
#		for j in range(0,Ly,1):
#			for k in range(0,Lz,1):
#				print i*j*k/(Lx*Ly*Lz)
		nlid = []
#		if (xyz[i,j,k] == 1):
		nlid+=[xyz_num[i,j,k]-1]
		n = int(xyz_num[i,j,k])
		mask = []
		for l in alist[n-1]:
			mask += [[l]]
		i1 = i-1
		i2 = i+1
		j1 = j-1
		j2 = j+1
		k1 = k-1
		k2 = k+1
		if (i1>=0):
			if (xyz[i1,j,k] == 1):
				mask[12-1].append(alist[xyz_num[i1,j,k]-1][13-1])
				mask[11-1].append(alist[xyz_num[i1,j,k]-1][14-1])
				mask[25-1].append(alist[xyz_num[i1,j,k]-1][28-1])
				mask[26-1].append(alist[xyz_num[i1,j,k]-1][27-1])
				mask[6-1].append(alist[xyz_num[i1,j,k]-1][9-1])
				mask[7-1].append(alist[xyz_num[i1,j,k]-1][8-1])
				mask[18-1].append(alist[xyz_num[i1,j,k]-1][21-1])
				mask[19-1].append(alist[xyz_num[i1,j,k]-1][20-1])
				mask[1-1].append(alist[xyz_num[i1,j,k]-1][4-1])
				mask[2-1].append(alist[xyz_num[i1,j,k]-1][3-1])
				nlid+=[xyz_num[i1,j,k]-1]
		if (i2<np.size(xyz,axis=0)):
			if (xyz[i2,j,k] == 1):
				mask[13-1].append(alist[xyz_num[i2,j,k]-1][12-1])
				mask[14-1].append(alist[xyz_num[i2,j,k]-1][11-1])
				mask[28-1].append(alist[xyz_num[i2,j,k]-1][25-1])
				mask[27-1].append(alist[xyz_num[i2,j,k]-1][26-1])
				mask[9-1].append(alist[xyz_num[i2,j,k]-1][6-1])
				mask[8-1].append(alist[xyz_num[i2,j,k]-1][7-1])
				mask[21-1].append(alist[xyz_num[i2,j,k]-1][18-1])
				mask[20-1].append(alist[xyz_num[i2,j,k]-1][19-1])
				mask[4-1].append(alist[xyz_num[i2,j,k]-1][1-1])
				mask[3-1].append(alist[xyz_num[i2,j,k]-1][2-1])
				nlid+=[xyz_num[i2,j,k]-1]
		if (j1>=0):
			if (xyz[i,j1,k] == 1):
				mask[1-1].append(alist[xyz_num[i,j1,k]-1][11-1])
				mask[2-1].append(alist[xyz_num[i,j1,k]-1][12-1])
				mask[16-1].append(alist[xyz_num[i,j1,k]-1][30-1])
				mask[5-1].append(alist[xyz_num[i,j1,k]-1][15-1])
				mask[17-1].append(alist[xyz_num[i,j1,k]-1][31-1])
				mask[4-1].append(alist[xyz_num[i,j1,k]-1][14-1])
				mask[3-1].append(alist[xyz_num[i,j1,k]-1][13-1])
				nlid+=[xyz_num[i,j1,k]-1]
		if (j2<np.size(xyz,axis=1)):
			if (xyz[i,j2,k] == 1):
				mask[11-1].append(alist[xyz_num[i,j2,k]-1][1-1])
				mask[12-1].append(alist[xyz_num[i,j2,k]-1][2-1])
				mask[30-1].append(alist[xyz_num[i,j2,k]-1][16-1])
				mask[15-1].append(alist[xyz_num[i,j2,k]-1][5-1])
				mask[31-1].append(alist[xyz_num[i,j2,k]-1][17-1])
				mask[14-1].append(alist[xyz_num[i,j2,k]-1][4-1])
				mask[13-1].append(alist[xyz_num[i,j2,k]-1][3-1])
				nlid+=[xyz_num[i,j2,k]-1]
		if (k1>=0):
			if (xyz[i,j,k1] == 1):
				mask[11-1].append(alist[xyz_num[i,j,k1]-1][12-1])
				mask[14-1].append(alist[xyz_num[i,j,k1]-1][13-1])
				mask[25-1].append(alist[xyz_num[i,j,k1]-1][26-1])
				mask[28-1].append(alist[xyz_num[i,j,k1]-1][27-1])
				mask[6-1].append(alist[xyz_num[i,j,k1]-1][7-1])
				mask[9-1].append(alist[xyz_num[i,j,k1]-1][8-1])
				mask[18-1].append(alist[xyz_num[i,j,k1]-1][19-1])
				mask[21-1].append(alist[xyz_num[i,j,k1]-1][20-1])
				mask[1-1].append(alist[xyz_num[i,j,k1]-1][2-1])
				mask[4-1].append(alist[xyz_num[i,j,k1]-1][3-1])
				nlid+=[xyz_num[i,j,k1]-1]
		if (k2<np.size(xyz,axis=2)):
			if (xyz[i,j,k2] == 1):
				mask[12-1].append(alist[xyz_num[i,j,k2]-1][11-1])
				mask[13-1].append(alist[xyz_num[i,j,k2]-1][14-1])
				mask[26-1].append(alist[xyz_num[i,j,k2]-1][25-1])
				mask[27-1].append(alist[xyz_num[i,j,k2]-1][28-1])
				mask[7-1].append(alist[xyz_num[i,j,k2]-1][6-1])
				mask[8-1].append(alist[xyz_num[i,j,k2]-1][9-1])
				mask[19-1].append(alist[xyz_num[i,j,k2]-1][18-1])
				mask[20-1].append(alist[xyz_num[i,j,k2]-1][21-1])
				mask[2-1].append(alist[xyz_num[i,j,k2]-1][1-1])
				mask[3-1].append(alist[xyz_num[i,j,k2]-1][4-1])
				nlid+=[xyz_num[i,j,k2]-1]
		if (i1>=0 and j1>=0):
			if (xyz[i1,j1,k] == 1):
				mask[1-1].append(alist[xyz_num[i1,j1,k]-1][14-1])
				mask[2-1].append(alist[xyz_num[i1,j1,k]-1][13-1])
				nlid+=[xyz_num[i1,j1,k]-1]
		if (i2<np.size(xyz,axis=0) and j1>=0):
			if (xyz[i2,j1,k] == 1):
				mask[3-1].append(alist[xyz_num[i2,j1,k]-1][12-1])
				mask[4-1].append(alist[xyz_num[i2,j1,k]-1][11-1])
				nlid+=[xyz_num[i2,j1,k]-1]
		if (i1>=0 and j2<np.size(xyz,axis=1)):
			if (xyz[i1,j2,k] == 1):
				mask[11-1].append(alist[xyz_num[i1,j2,k]-1][4-1])
				mask[12-1].append(alist[xyz_num[i1,j2,k]-1][3-1])
				nlid+=[xyz_num[i1,j2,k]-1]
		if (i2<np.size(xyz,axis=0) and j2<np.size(xyz,axis=1)):
			if (xyz[i2,j2,k] == 1):
				mask[14-1].append(alist[xyz_num[i2,j2,k]-1][1-1])
				mask[13-1].append(alist[xyz_num[i2,j2,k]-1][2-1])
				nlid+=[xyz_num[i2,j2,k]-1]
		if (i1>=0 and k1>=0):
			if (xyz[i1,j,k1] == 1):
				mask[11-1].append(alist[xyz_num[i1,j,k1]-1][13-1])
				mask[25-1].append(alist[xyz_num[i1,j,k1]-1][27-1])
				mask[6-1].append(alist[xyz_num[i1,j,k1]-1][8-1])
				mask[18-1].append(alist[xyz_num[i1,j,k1]-1][20-1])
				mask[1-1].append(alist[xyz_num[i1,j,k1]-1][3-1])
				nlid+=[xyz_num[i1,j,k1]-1]
		if (i2<np.size(xyz,axis=0) and k1>=0):
			if (xyz[i2,j,k1] == 1):
				mask[14-1].append(alist[xyz_num[i2,j,k1]-1][12-1])
				mask[28-1].append(alist[xyz_num[i2,j,k1]-1][26-1])
				mask[9-1].append(alist[xyz_num[i2,j,k1]-1][7-1])
				mask[21-1].append(alist[xyz_num[i2,j,k1]-1][19-1])
				mask[4-1].append(alist[xyz_num[i2,j,k1]-1][2-1])
				nlid+=[xyz_num[i2,j,k1]-1]
		if (i1>=0 and k2<np.size(xyz,axis=2)):
			if (xyz[i1,j,k2] == 1):
				mask[12-1].append(alist[xyz_num[i1,j,k2]-1][14-1])
				mask[26-1].append(alist[xyz_num[i1,j,k2]-1][28-1])
				mask[7-1].append(alist[xyz_num[i1,j,k2]-1][9-1])
				mask[19-1].append(alist[xyz_num[i1,j,k2]-1][21-1])
				mask[2-1].append(alist[xyz_num[i1,j,k2]-1][4-1])
				nlid+=[xyz_num[i1,j,k2]-1]
		if (i2<np.size(xyz,axis=0) and k2<np.size(xyz,axis=2)):
			if (xyz[i2,j,k2] == 1):
				mask[13-1].append(alist[xyz_num[i2,j,k2]-1][11-1])
				mask[27-1].append(alist[xyz_num[i2,j,k2]-1][25-1])
				mask[8-1].append(alist[xyz_num[i2,j,k2]-1][6-1])
				mask[20-1].append(alist[xyz_num[i2,j,k2]-1][18-1])
				mask[3-1].append(alist[xyz_num[i2,j,k2]-1][1-1])
				nlid+=[xyz_num[i2,j,k2]-1]
		if (j1>=0 and k1>=0):
			if (xyz[i,j1,k1] == 1):
				mask[1-1].append(alist[xyz_num[i,j1,k1]-1][12-1])
				mask[4-1].append(alist[xyz_num[i,j1,k1]-1][13-1])
				nlid+=[xyz_num[i,j1,k1]-1]
		if (j2<np.size(xyz,axis=1) and k1>=0):
			if (xyz[i,j2,k1] == 1):
				mask[11-1].append(alist[xyz_num[i,j2,k1]-1][2-1])
				mask[14-1].append(alist[xyz_num[i,j2,k1]-1][3-1])
				nlid+=[xyz_num[i,j2,k1]-1]
		if (j1>=0 and k2<np.size(xyz,axis=2)):
			if (xyz[i,j1,k2] == 1):
				mask[2-1].append(alist[xyz_num[i,j1,k2]-1][11-1])
				mask[3-1].append(alist[xyz_num[i,j1,k2]-1][14-1])
				nlid+=[xyz_num[i,j1,k2]-1]
		if (j2<np.size(xyz,axis=1) and k2<np.size(xyz,axis=2)):
			if (xyz[i,j2,k2] == 1):
				mask[12-1].append(alist[xyz_num[i,j2,k2]-1][1-1])
				mask[13-1].append(alist[xyz_num[i,j2,k2]-1][4-1])
				nlid+=[xyz_num[i,j2,k2]-1]
		if (i1>=0 and j1>=0 and k1>=0):
			if (xyz[i1,j1,k1] == 1):
				mask[1-1].append(alist[xyz_num[i1,j1,k1]-1][13-1])
				nlid+=[xyz_num[i1,j1,k1]-1]
		if (i2<np.size(xyz,axis=0) and j1>=0 and k1>=0):
			if (xyz[i2,j1,k1] == 1):
				mask[4-1].append(alist[xyz_num[i2,j1,k1]-1][12-1])
				nlid+=[xyz_num[i2,j1,k1]-1]
		if (i1>=0 and j2<np.size(xyz,axis=1) and k1>=0):
			if (xyz[i1,j2,k1] == 1):
				mask[11-1].append(alist[xyz_num[i1,j2,k1]-1][3-1])
				nlid+=[xyz_num[i1,j2,k1]-1]
		if (i1>=0 and j1>=0 and k2<np.size(xyz,axis=2)):
			if (xyz[i1,j1,k2] == 1):
				mask[2-1].append(alist[xyz_num[i1,j1,k2]-1][14-1])
				nlid+=[xyz_num[i1,j1,k2]-1]
		if (i2<np.size(xyz,axis=0) and j2<np.size(xyz,axis=1) and k1>=0):
			if (xyz[i2,j2,k1] == 1):
				mask[14-1].append(alist[xyz_num[i2,j2,k1]-1][2-1])
				nlid+=[xyz_num[i2,j2,k1]-1]
		if (i1>=0 and j2<np.size(xyz,axis=1) and k2<np.size(xyz,axis=2)):
			if (xyz[i1,j2,k2] == 1):
				mask[12-1].append(alist[xyz_num[i1,j2,k2]-1][4-1])
				nlid+=[xyz_num[i1,j2,k2]-1]
		if (i2<np.size(xyz,axis=0) and j1>=0 and k2<np.size(xyz,axis=2)):
			if (xyz[i2,j1,k2] == 1):
				mask[3-1].append(alist[xyz_num[i2,j1,k2]-1][11-1])
				nlid+=[xyz_num[i2,j1,k2]-1]
		if (i2<np.size(xyz,axis=0) and j2<np.size(xyz,axis=1) and k2<np.size(xyz,axis=2)):
			if (xyz[i2,j2,k2] == 1):
				mask[13-1].append(alist[xyz_num[i2,j2,k2]-1][1-1])
				nlid+=[xyz_num[i2,j2,k2]-1]
		nlid2 = []
		for q in nlid:
			nlid2 = np.concatenate((nlid2,alist2[q]))
		for q in range(len(mask)):
			for r in range(len(mask[q])):
				for m in nlid2:
					for o in range(len(nlist[int(m)-1])):
						for p in range(len(nlist[int(m)-1][o])):
							if (nlist[int(m)-1][o][p]==mask[q][r]):
								nlist[int(m)-1][o][p] = min(mask[q])
				for m in nlid:
					for o in range(len(alist[int(m)])):
						if (alist[int(m)][o]==mask[q][r]):
							alist[int(m)][o] = min(mask[q])
	return alist,nlist


# In[17]:


def dist(xyz,xyz_num,nuc,ifxyz,ifxyz_num,ifnuc,gxyz,gxyz_num,gnuc,La,Lb,Lc):
	Lpbca = np.size(xyz,axis=0)*La
	Lpbcb = np.size(xyz,axis=1)*Lb
	Lpbcc = np.size(xyz,axis=2)*Lc
	print Lpbca,Lpbcb,Lpbcc
	d_alist = np.zeros((3,31*nuc+15*ifnuc+15*gnuc))
	for i in np.ndindex(xyz.shape):
		if xyz[i]==1:
			j = 31*(xyz_num[i]-1)
			d_alist[0,j]=d_alist[0,j+17]=d_alist[0,j+5]=d_alist[0,j+24]=d_alist[0,j+10]=0+La*i[0]
			d_alist[0,j+4]=d_alist[0,j+21]=d_alist[0,j+9]=d_alist[0,j+28]=d_alist[0,j+14]=La*0.5+La*i[0]
			d_alist[0,j+1]=d_alist[0,j+18]=d_alist[0,j+6]=d_alist[0,j+25]=d_alist[0,j+11]=0+La*i[0]
			d_alist[0,j+3]=d_alist[0,j+20]=d_alist[0,j+8]=d_alist[0,j+27]=d_alist[0,j+13]=La+La*i[0]
			d_alist[0,j+2]=d_alist[0,j+19]=d_alist[0,j+7]=d_alist[0,j+26]=d_alist[0,j+12]=La+La*i[0]
			d_alist[0,j+15]=d_alist[0,j+22]=d_alist[0,j+29]=La*0.25+La*i[0]
			d_alist[0,j+16]=d_alist[0,j+23]=d_alist[0,j+30]=La*0.75+La*i[0]

			d_alist[1,j:j+5] = d_alist[1,j+15] = d_alist[1,j+16] = 0+Lb*i[1]
			d_alist[1,j+5:j+10] = d_alist[1,j+22] = d_alist[1,j+23] = Lb*0.5+Lb*i[1]
			d_alist[1,j+10:j+15] = d_alist[1,j+29] = d_alist[1,j+30] = Lb+Lb*i[1]
			d_alist[1,j+17:j+22] = Lb*0.25+Lb*i[1]
			d_alist[1,j+24:j+29] = Lb*0.75+Lb*i[1]

			d_alist[2,j]=d_alist[2,j+17]=d_alist[2,j+5]=d_alist[2,j+24]=d_alist[2,j+10]=d_alist[2,j+3]=d_alist[2,j+20]=d_alist[2,j+8]=d_alist[2,j+27]=d_alist[2,j+13]=0+Lc*i[2]
			d_alist[2,j+1]=d_alist[2,j+18]=d_alist[2,j+6]=d_alist[2,j+25]=d_alist[2,j+11]=d_alist[2,j+2]=d_alist[2,j+19]=d_alist[2,j+7]=d_alist[2,j+26]=d_alist[2,j+12]=Lc+Lc*i[2]
			d_alist[2,j+15]=d_alist[2,j+16]=d_alist[2,j+29]=d_alist[2,j+30]=0.25*Lc+Lc*i[2]
			d_alist[2,j+22]=d_alist[2,j+23]=0.75*Lc+Lc*i[2]
			d_alist[2,j+4]=d_alist[2,j+21]=d_alist[2,j+9]=d_alist[2,j+28]=d_alist[2,j+14]=0.5*Lc+Lc*i[2]

	for i in np.ndindex(ifxyz.shape):
		if ifxyz[i]==1:
			j = 15*(ifxyz_num[i]-1 - nuc) + 31*nuc
			d_alist[0,j]=d_alist[0,j+5]=d_alist[0,j+10]=0+La*i[0]
			d_alist[0,j+4]=d_alist[0,j+9]=d_alist[0,j+14]=La*0.5+La*i[0]
			d_alist[0,j+1]=d_alist[0,j+6]=d_alist[0,j+11]=0+La*i[0]
			d_alist[0,j+3]=d_alist[0,j+8]=d_alist[0,j+13]=La+La*i[0]
			d_alist[0,j+2]=d_alist[0,j+7]=d_alist[0,j+12]=La+La*i[0]

			d_alist[1,j:j+5] = 0+Lb*i[1]
			d_alist[1,j+5:j+10] = Lb*0.5+Lb*i[1]
			d_alist[1,j+10:j+15] = Lb+Lb*i[1]

			d_alist[2,j]=d_alist[2,j+5]=d_alist[2,j+10]=d_alist[2,j+3]=d_alist[2,j+8]=d_alist[2,j+13]=0+Lc*i[2]
			d_alist[2,j+1]=d_alist[2,j+6]=d_alist[2,j+11]=d_alist[2,j+2]=d_alist[2,j+7]=d_alist[2,j+12]=Lc+Lc*i[2]
			d_alist[2,j+4]=d_alist[2,j+9]=d_alist[2,j+14]=0.5*Lc+Lc*i[2]

	for i in np.ndindex(gxyz.shape):
		if gxyz[i]==1:
			j = 15*(gxyz_num[i]-1 - nuc - ifnuc) + 31*nuc + 15*ifnuc
			d_alist[0,j]=d_alist[0,j+5]=d_alist[0,j+10]=0+La*i[0]
			d_alist[0,j+4]=d_alist[0,j+9]=d_alist[0,j+14]=La*0.5+La*i[0]
			d_alist[0,j+1]=d_alist[0,j+6]=d_alist[0,j+11]=0+La*i[0]
			d_alist[0,j+3]=d_alist[0,j+8]=d_alist[0,j+13]=La+La*i[0]
			d_alist[0,j+2]=d_alist[0,j+7]=d_alist[0,j+12]=La+La*i[0]

			d_alist[1,j:j+5] = 0+Lb*i[1]
			d_alist[1,j+5:j+10] = Lb*0.5+Lb*i[1]
			d_alist[1,j+10:j+15] = Lb+Lb*i[1]

			d_alist[2,j]=d_alist[2,j+5]=d_alist[2,j+10]=d_alist[2,j+3]=d_alist[2,j+8]=d_alist[2,j+13]=0+Lc*i[2]
			d_alist[2,j+1]=d_alist[2,j+6]=d_alist[2,j+11]=d_alist[2,j+2]=d_alist[2,j+7]=d_alist[2,j+12]=Lc+Lc*i[2]
			d_alist[2,j+4]=d_alist[2,j+9]=d_alist[2,j+14]=0.5*Lc+Lc*i[2]
	np.save(sys.argv[1]+'_coordinates',d_alist)
	del d_alist


# In[18]:


def merge_neighbors():
	global neighbors,num_neighbors
	j = -1
	for i in neighbors[:,0]:
		j = j+1
		index = (neighbors[j]!=0)[1:]
		index = np.append([False],index)
		neighbors[j,index] = neighbors[neighbors[j,index]-1,0]
		del index
	print "Merged neighborlist"
	rows = np.unique(neighbors[:,0])
	unique_neighbors = np.memmap(sys.argv[1]+"_neighbors.npy",dtype=int,mode='w+',shape=(len(rows),7))
	np.save(sys.argv[1]+"_shape",[len(rows),7])
	#unique_neighbors[:,:] = -1
	unique_neighbors[:,0] = rows[:]
	ID = {}
	j = -1
	for i in rows:
		j=j+1
		ID[i]=j
	del rows
	j = -1
	print "Unique empty neighborlist created"
	for i in neighbors[:,0]:
		j = j+1
		k = ID[i]#np.where(unique_neighbors[:,0]==i)[0][0]
		arr1 = neighbors[j,1:]
		arr2 = unique_neighbors[k,1:]
		arr = np.union1d(arr1[arr1!=0],arr2[arr2!=0])
		arr = np.append(arr,[0]*(len(arr2)-len(arr)))
		unique_neighbors[k,1:len(arr)+1] = arr[:]
	print "Unique neighborlist created"
	del neighbors
	os.remove('testmmp_2.npy')
	print "Old neighborlist removed"
	del num_neighbors
	num_neighbors = 0*unique_neighbors[:,0]
	print "Number of neighbors"
	j = -1
	for i in unique_neighbors:
		j = j+1
		num_neighbors[j] = len(unique_neighbors[j,unique_neighbors[j]!=0])-1
	global neighbors
	neighbors = unique_neighbors
	del unique_neighbors
	print np.unique(num_neighbors)
	print neighbors[neighbors[:,0]==12]
	print neighbors
	#for i in neighbors[:,0]:
	#	j = j+1
	#	arr1 = 
	#	arr = np.union1d(neighbors[j,1:len(neighbors[j,neighbors[j]!=-1])],neighbors[i-1,1:len(neighbors[i-1,neighbors[i-1]!=-1])])
	#	neighbors[i-1,1:len(arr)+1] = arr[:]

	#rows = np.unique(neighbors[:,0])[:]-1
	#unique_neighbors = np.memmap('testmmp_unique.npy',dtype=int,mode='w+',shape=(len(rows),20))
	#unique_neighbors[:,:20] = neighbors[rows[:],:20]
	#del neighbors
	#os.remove('testmmp.npy')
	#print unique_neighbors
	#del num_neighbors
	#num_neighbors = 0*unique_neighbors[:,0]
	#j = -1
	#for i in unique_neighbors:
	#	j = j+1
	#	num_neighbors[j] = len(unique_neighbors[j,unique_neighbors[j]!=-1])-1
	#	if (num_neighbors[j]>6):
	#		print unique_neighbors[j]
	#print num_neighbors
	#print np.unique(num_neighbors)
	#neighbors = np.memmap('testmmp.npy',dtype=int,mode='r',shape=(rows,10))


# In[19]:


def nonint_neighborlist(start,end):
	global siteID,neighbors
	
	for i in range(len(neighbors)):
		if (len(neighbors[i])>1):
			for j in neighbors[i][1]:
				if (len(neighbors[j-1])>1):
					neighbors[j-1][1]+=[neighbors[i][0][0]]
				if (len(neighbors[j-1])==1):
					neighbors[j-1].append([neighbors[i][0][0]])
	#return alist,nlist#nlist2


# In[20]:


def mapmerge(i):
	global nlist3,nlist2
	for j in range(len(nlist)):
		if (nlist3[i][0][0]==nlist[j][0][0]):
			nlist3[i][1]=list(set().union(nlist3[i][1],nlist[j][1]))
	return nlist3[i]


# In[21]:


def mergelist(alist,nlist):
	#global nlist3
	#global nlist2
	#nlist2 = copy.deepcopy(nlist)
	nlist3 = []
	ualist = []
	for i in alist:
		for j in i:
			ualist += [j]
#	print "Alist copy created"
	ualist = np.unique(np.asarray(ualist))
#	print "Unique alist created"
	rows = np.size(ualist)
	ID = {}
	for row in range(rows):
		nlist3+=[[[ualist[row]*1],[]]]
		ID[ualist[row]]=row
#	print "Unique empty nlist created"
	#p = Pool(480)
	#result = p.map(mapmerge,range(len(nlist3)))
	#for i in range(len(nlist3)):
	#	for j in range(len(nlist)):
	#		if (nlist3[i][0][0]==nlist[j][0][0]):
	#			nlist3[i][1]=list(set().union(nlist3[i][1],nlist[j][1]))
	for i in range(len(nlist)):
		j = ID[nlist[i][0][0]]
		nlist3[j][1]=list(set().union(nlist3[j][1],nlist[i][1]))
	print "Merged"
	return alist,nlist3


# In[22]:


def animate(i,simstep,dlist,points):
	for j,k in zip(simstep[:,i],points):
		data = dlist[:,j-1:j]
		k.set_data(data[0:2])
		k.set_3d_properties(data[2])
	return points


# In[23]:


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


# In[24]:


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


# In[25]:


def merging():
	global siteID,neighbors
	nlist3 = []
	ualist = []
	for i in siteID:
		for j in i:
			ualist += [j]
	ualist = np.unique(np.asarray(ualist))
	rows = np.size(ualist)
	ID = {}
	row = 0
	while row < rows:
		nlist3+=[[[ualist[row]*1],[]]]
		ID[ualist[row]]=row
		row = row+1
	del ualist
	i = 0
	while i < len(neighbors):
		j = ID[neighbors[i][0][0]]
		nlist3[j][1]=list(set().union(nlist3[j][1],neighbors[i][1]))
		i = i+1
	print "Merged"
	del neighbors
	return nlist3


# ## MAIN ####

# In[ ]:

# Import the generated unit cell map
# _z : Unit cell of zeolite in the system
# _ig : Unit cell of the interface between zeolite and external gas, composed of a single layer of unit cells
# _g : Unit cell of the gas phase, obtained as difference between system (all unit cells) - zeolite - interface

lfile_z = sys.argv[1]+"_z.npy"
lfile_ig = sys.argv[1]+"_ig.npy"
lfile_g = sys.argv[1]+"_g.npy"

# Size of the unit cell in x = La, y = Lb, z = Lc

La = 2.01
Lb = 1.99
Lc = 1.34

# Global parameters that contain information about sites, neighbors of sites, and each type of site: Zeolite type, Interface type, Gas phase type

siteID = []
neighbors = []
Zeo=[]
Zeo_num=[]
ZG=[]
ZG_num=[]
Gas=[]
Gas_num=[]

# Pseudo length to avoid jump across sites separated by a periodic boundary (if present)
# The factor 1.2 is greater than the actual distance between sites, i.e. factor of 1.0 and smaller than the next equivalent site, i.e. factor of 2.0
# This ensures that only jumps that cross a periodic boundary are wrapped, while if the coordinates are allowed to be at periodic boundary (if present)

Llist = {}
Llist['ZI'] = {}
Llist['ZZ'] = {}
Llist['ZS'] = {}
Llist['ZGI']= {}
Llist['GI'] = {}
Llist['ZI']['ZZ'] = {}
Llist['ZI']['ZZ']['x'] = 1.2*0.5*La
Llist['ZI']['ZZ']['y'] = 1.2*0.5*Lb
Llist['ZI']['ZZ']['z'] = 1.2*0.5*Lc
Llist['ZZ'] = {}
Llist['ZZ']['ZI']   = {}
Llist['ZZ']['ZI']['x'] = 1.2*0.5*La
Llist['ZZ']['ZI']['y'] = 1.2*0.5*Lb
Llist['ZZ']['ZI']['z'] = 1.2*0.5*Lc

Llist['ZI']['ZS']        = {}
Llist['ZI']['ZS']['x']   = 1.2*0.5*La
Llist['ZI']['ZS']['y']   = 1.2*0.5*Lb
Llist['ZI']['ZS']['z']   = 1.2*0.5*Lc
Llist['ZS']['ZI']        = {}
Llist['ZS']['ZI']['x']   = 1.2*0.5*La
Llist['ZS']['ZI']['y']   = 1.2*0.5*Lb
Llist['ZS']['ZI']['z']   = 1.2*0.5*Lc

Llist['ZI']['ZGI']       = {}
Llist['ZI']['ZGI']['x']  = 1.2*La
Llist['ZI']['ZGI']['y']  = 1.2*Lb
Llist['ZI']['ZGI']['z']  = 1.2*Lc

Llist['ZGI']['ZI']       = {}
Llist['ZGI']['ZI']['x']  = 1.2*La
Llist['ZGI']['ZI']['y']  = 1.2*Lb
Llist['ZGI']['ZI']['z']  = 1.2*Lc

Llist['ZGI']['ZGI'] = {}
Llist['ZGI']['ZGI']['x'] = 1.2*La
Llist['ZGI']['ZGI']['y'] = 1.2*Lb
Llist['ZGI']['ZGI']['z'] = 1.2*Lc

Llist['ZGI']['GI'] = {}
Llist['ZGI']['GI']['x']  = 1.2*La
Llist['ZGI']['GI']['y']  = 1.2*Lb
Llist['ZGI']['GI']['z']  = 1.2*Lc

Llist['GI']['ZGI'] = {}
Llist['GI']['ZGI']['x']  = 1.2*La
Llist['GI']['ZGI']['y']  = 1.2*Lb
Llist['GI']['ZGI']['z']  = 1.2*Lc

Llist['GI']['GI'] = {}
Llist['GI']['GI']['x']   = 1.2*La
Llist['GI']['GI']['y']   = 1.2*Lb
Llist['GI']['GI']['z']   = 1.2*Lc
np.save(sys.argv[1]+"_length",Llist)
del Llist

# Initialize the generation of non-periodic mfi-type unit cell sites where zeolite type unit cell is present
# Initialize the generation of non-periodic, only intersection sites, of the mfi type unit cell where gas and interface are present
# Gas and interface are not restricted by porous constraints of zeolite, and therefore are initialized by have all possible connections with the neighbor intersection sites. 
# --> an adsorbing molecule can move from interface to zeolite (and vice-versa) through paths between intersection type sites of the two unit cells

initialize_z(lfile_z,lfile_ig,lfile_g)#,siteID,neighbors,label)
print "\n Initialized \n"

# Remove the connections in the interface unit cell according to the neighbor zeolite unit cell to ensure continuity of connections.  Note, each unit cell has 26 neighbors, hence 26 conditions

remove_neighborlist_if_z()

# Add some connections to ensure conitnuity of connections with the zeolite lattice.  Note, each unit cell has 26 neighbors, hence 26 conditions

add_neighborlist_if_z()

# Merge the unit cells together by identifying which sites are the same across neighboring unit cells. Note, each unit cell has 26 neighbors, hence 26 conditions

nonperiodic_neighborlist_2()
print len(np.unique(neighbors[:,0]))
print int((time.time()-t0)/60),'min'
print neighbors

# Merge the unit cells of the gas-phase to include periodicity. This ability was present but unused in the simulations

periodic_neighborlist_g()
print neighbors

# Merge the duplicate sites of nerighboring unit cells and combine the neighborlist that represent connected paths where the molecule can hop to from that site

merge_neighbors()
print np.unique(num_neighbors)
del neighbors

# Create a coordinate list of the sites of this lattice according to the unit cell sizes of specfied above with La, Lb, Lc

dist(Zeo,Zeo_num,nZuc,ZG,ZG_num,nZGuc,Gas,Gas_num,nGuc,La,Lb,Lc)

