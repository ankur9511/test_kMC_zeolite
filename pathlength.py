#!/usr/bin/env python2
# coding: utf-8

# In[5]:


import numpy as np
import glob
from matplotlib import rcParams
params = {'backend': 'agg'}
rcParams['font.family'] = 'serif'
rcParams.update(params)
import matplotlib.style
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cfit
from scipy import integrate
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import interp1d
import scipy.stats
import sys
from scipy.signal import savgol_filter


# In[ ]:


itr = 'test'
nxwindow = 10
nxoverlap = 10
fin_name = ["."]
fin_label = [r'$Test$']
marker = ["s"]
hloc = ['left']
zorder = [0]
cf = ['#1F77B4']
norm = [3]
yloc_r_ax1 = -5.0
yloc_r_ax2 = -5.0
xmin_ax1 = 1
xmax_ax1 = 10**7 #5*10**4
xmin_ax2 = 1
xmax_ax2 = 10**7 #5*10**4
xleg_ax1 = 3*10**3
xleg_ax2 = 3*10**3
ymin_ax2 = -18
ymin_ax1 = -18
ymax_ax1 = -2
ymax_ax2 = -2
flag = 1


# ### Input parameters for the plot
# ##### opt: Flag to different between many plotting styles, if needed
# ##### itr: A string used for defaut naming of output from this iteration of code
# ##### nxwindow: 
# ##### nxoverlap:
# ##### fin_name: Location of input data in specified order
# ##### fin_label: Label for the line plot based on input data order
# ##### marker: Type of data point marker
# ##### zorder: Order of line on top of another
# ##### cf: Color of marker and line plot

# ## Functions of the code

# ### Get names of all trajectory files at the location

# In[9]:


def getfiles(dirname):
	names = glob.glob(dirname+"/traj*_*")
	return names


# ### Calculate path length for a trajectory file

# In[12]:


def readfile(name,flag):
	f = open(name,'r')
	i = 0
	flag = 0
	pathlength = np.nan
	lines = reversed(f.readlines())
	for line in lines:
		i = i+1
		if (i==1 and line=="Finished"):
			flag = 1
		secondlastline = "0"
		if (i==3 and flag==1):
			secondlastline = line
			if (len(secondlastline.split("\t"))==5):
				pathlength = secondlastline.split("\t")[4]
	#if float(pathlength) < 0.995:
	#	pathlength = np.nan
	f.close()
	return float(pathlength)


# ### Wrap path length calculation for every trajectory file

# In[11]:


def escape_time_cumuprob(dirname,flag):
	names = getfiles(dirname)
	t = []
	for name in names:
		t+=[readfile(name,flag)]
	t = np.array(t)
	t = t[np.isfinite(t)]
	return t


# ### Create histogram

# In[10]:


def overlaphistogram(xy,xwindow,xoverlap,dx):
	xedge = np.exp(np.arange(np.log(0.1),np.log(xy.max())+1.2,0.4))
	#xedge2 = np.exp(np.arange(np.log(xedge1[-1]),np.log(xy.max())+0.25,0.25))
	#xedge = np.concatenate((xedge1[:-1],xedge2))
	ynew,xedge = np.histogram(xy,bins=xedge)
	xcenter = (xedge[1:]+xedge[:-1])*0.5
	del xedge
	return xcenter,ynew


# ### Plot histogram

# In[5]:


def plotcumtime(ax1,scale,hloc,ni,dirname,lab,color,marker,zorder):
	global itr,nxoverlap,nxwindow,fin_name,fin_label,nofin_name,nofin_label,flag,lines,labels,yloc_r_ax2
	xy = escape_time_cumuprob(dirname,flag)
	xy = xy[:min(144000,len(xy))]

	###### Plot 1: Cumulative probability ####
	x,y = overlaphistogram(xy,nxwindow,nxoverlap,0.1)
	avgt = np.average(xy)
	avgterr = np.sqrt(np.var(xy)/len(xy))
	y_dash = y/integrate.simps(y,x)
	y_dash = np.log(y_dash)
	print "y-values are: ",y_dash
	print "x-values are: ",x
	l1,=ax1.semilogx(x,y_dash,linewidth=0,label=lab,color=color,marker=marker,markersize=2,zorder=zorder)
	ax1.semilogx(x,y_dash,linewidth=2,markersize=0,alpha=0.5,color=color,zorder=7)
	ax1.semilogx([avgt,avgt],[ymin_ax2,ymax_ax2],color=color,linestyle='--',linewidth=0.5)
	ax1.text(avgt,yloc_r_ax2,r'$\mathbf{\/\/\/\/\/\/\/\/\langle r \rangle\/=\/}$'+"{:3.0f}      ".format(int(avgt)),fontsize=6,horizontalalignment=hloc,weight='bold',color=color)
	#yloc_r_ax2 += -0.25
	ax1.set_xlabel(r'$\mathrm{r\/(nm)}$')
	ax1.set_ylabel(r'$\mathrm{ln[P(r)]}$')
	lines[0]+=[l1]
	labels[0]+=[lab]
	return x,y_dash


# ## Main part of the Code

# In[ ]:


plt.rc('font', family='serif')
fig1 = plt.figure(figsize = (3.5,3.5))
ax1 = fig1.add_subplot(111)
lines = [[]]
labels = [[]]
miny = 10
maxx = 0
#print len(fin_name),len(axlist),len(hloc),len(norm),len(fin_label),len(cf),len(marker),len(zorder)
for i in range(len(fin_name)):
	x,y = plotcumtime(ax1,1.0,hloc[i],norm[i],dirname=fin_name[i],lab=fin_label[i],color=cf[i],marker=marker[i],zorder=zorder[i])
	maxx = max(x[-1],maxx)
	miny = min(np.amin(y[np.isfinite(y)]),miny)

ax1.set_xlim(xmin=0.8,xmax=maxx)
ax1.set_ylim(ymin=ymin_ax2,ymax=ymax_ax2)
ax1.set_xticks([10**1,10**3,10**5])#,10**-2,10**0,10**2,10**4]

leg1 = fig1.legend(lines[0],labels[0],loc=[0.65,0.65],fancybox=False,framealpha=None,frameon=False,fontsize=6)

for line, text in zip(leg1.get_lines(), leg1.get_texts()):
    text.set_color(line.get_color())

# fig1.text(.5,-0.1,'Fig: Path lengths obtained for n_sim = 70000',ha='center') 

fig1.savefig("Path_lengths"+itr+".png",dpi=1000,bbox_inches='tight')
fig1.savefig("Path_lengths"+itr+".eps",dpi=1000,bbox_inches='tight')
fig1.savefig("Path_lengths"+itr+".pdf",dpi=1000,bbox_inches='tight')

