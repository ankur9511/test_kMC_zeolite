#!/usr/bin/env python2
# coding: utf-8

# ##### python rel_surface_area.py 50_50_4_500_z.npy,30_300_v2/50_50_10_500_z.npy,30_300_v2/50_50_25_500_z.npy,30_300_v2/50_50_50_500_z.npy 4,10,25,50 0.00033,0.00026,0.00015,0.00010 0_0_0_500new_z.npy

# In[ ]:


import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import rc
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


arr = sys.argv[1].split(',')
labarr = sys.argv[2].split(',')
Narr = sys.argv[3].split(',')
a = np.load(sys.argv[4])
USAx = 1.99*1.34
USAy = 2.01*1.34
USAz = 1.99*2.01
SAx = surfarea0(a,USAx)
SAy = surfarea1(a,USAy)
SAz = surfarea2(a,USAz)
print "For system %s" % (sys.argv[3])
print "SAx = ",SAx," SAy = ",SAy," SAz = ",SAz," nm^2"
refTSA  = SAx+SAy+SAz
print "TSA = ",refTSA," nm^2"
#print "Reference area",6*300*300," nm^2"
refTV = np.sum(a)*2.01*1.99*1.34
print "Volume",refTV," nm^3"
refSpecA = refTSA/refTV
print "SSA in m^2/g : %.1f " % (((refTSA/refTV)/(10**(-9)))/(1787587.597))
fig = plt.figure(figsize=(3,3))
ax = fig.add_axes([0.05,0.05,0.95,0.95])
#ax.text(600,200,r'$\alpha/\gamma ,\/\alpha,\/\gamma\/N$',style='normal')
for aname,label,N in zip(arr,labarr,Narr):
	a = np.load(aname)
	USAx = 1.99*1.34
	USAy = 2.01*1.34
	USAz = 1.99*2.01
	SAx = surfarea0(a,USAx)
	SAy = surfarea1(a,USAy)
	SAz = surfarea2(a,USAz)
	#print "For system %s" % (aname.split('/')[-1])
	#print "SAx = ",SAx," SAy = ",SAy," SAz = ",SAz," nm^2"
	tTSA  = SAx+SAy+SAz
	#print "TSA = ",tTSA," nm^2"
	#print "Reference area",6*300*300," nm^2"
	tTV = np.sum(a)*2.01*1.99*1.34
	#print "Volume",tTV," nm^3"
	#print "Reference Volume",300*300*300," nm^3"
	print "SSA in m^2/g : %.1f " % (((tTSA/tTV)/(10**(-9)))/(1787587.597)) #tTSA/tTV is 1/nm = 10^9 1/m , * m^3/g
	tSpecA = ((tTSA/tTV)-refSpecA)*100.0/refSpecA
	print "Increased Specific External Area for %s = %.1f" % (label,tSpecA)

	ax.plot(500,tSpecA,label=label,marker='*',markersize=5)
	ax.text(600,tSpecA, r'$\alpha/\gamma = %s ,\/\alpha = %s nm,\/\gamma = %s nm, \/ N = %s nm^{-2} $' % (str(round(50.0/float(label),1)),'50',str(int(label)),N),style='normal',fontsize=6)
ax.set_xlim(0,2000)
ax.set_ylim(0,200)
ax.set_yticks([0,20,40,60,80,100,120,140,160,180])
ax.set_xticks([0,500,1000,1500,2000])
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
ax.set_xlabel(r'$\beta$ (nm)',fontsize=8)
ax.set_ylabel('Increased Specific External Area (%)',fontsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.spines['bottom'].set_visible(True)
#ax.spines['left'].set_visible(True)
plt.savefig('rel_surf_spec_area.png',dpi=800,bbox_inches='tight')
plt.savefig('rel_surf_spec_area.pdf',dpi=800,bbox_inches='tight')
plt.savefig('rel_surf_spec_area.eps',dpi=800,bbox_inches='tight')
plt.savefig('rel_surf_spec_area.svg',dpi=800,bbox_inches='tight')

plt.show()

