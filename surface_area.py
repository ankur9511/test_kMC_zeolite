#!/usr/bin/env python2
# coding: utf-8

# In[ ]:


import numpy as np
import sys


# In[ ]:


def surfarea0(m,USAx):
	SAx = 0.0
	for k in range(np.size(m,axis=2)):
		for j in range(np.size(m,axis=1)):
			flag = m[0,j,k]
			for i in range(np.size(m,axis=0)):
				if m[i,j,k] != flag :
					SAx += USAx
					flag = m[i,j,k]
	return SAx


# In[ ]:


def surfarea1(m,USAy):
	SAy = 0.0
	for k in range(np.size(m,axis=2)):
		for i in range(np.size(m,axis=0)):
			flag = m[i,0,k]
			for j in range(np.size(m,axis=1)):
				if m[i,j,k] != flag :
					SAy += USAy
					flag = m[i,j,k]
	return SAy


# In[ ]:


def surfarea2(m,USAz):
	SAz = 0.0
	for i in range(np.size(m,axis=0)):
		for j in range(np.size(m,axis=1)):
			flag = m[i,j,0]
			for k in range(np.size(m,axis=2)):
				if m[i,j,k] != flag :
					SAz += USAz
					flag = m[i,j,k]
	return SAz


# In[ ]:


a = np.load(sys.argv[1])
USAx = 1.99*1.34
USAy = 2.01*1.34
USAz = 1.99*2.01
SAx = surfarea0(a,USAx)
SAy = surfarea1(a,USAy)
SAz = surfarea2(a,USAz)
print "SAx = ",SAx," SAy = ",SAy," SAz = ",SAz," nm^2"
print "TSA = ",SAx+SAy+SAz," nm^2"
#print "Reference area",6*300*300," nm^2"


# In[ ]:


print "Volume",np.sum(a)*2.01*1.99*1.34," nm^3"
#print "Reference Volume",300*300*300," nm^3"

