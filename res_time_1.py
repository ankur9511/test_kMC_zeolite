#!/usr/bin/env python2
# coding: utf-8

# In[2]:


import numpy as np
import glob
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


opt = sys.argv[1]
if (opt == '1'):
	itr = 'surf_high_v7'
	nxwindow = 10
	nxoverlap = 10
	fin_name = ["0_300/any_any_high","30_300_v2/30_30_300/any_any_high","30_300_v2/30_10_300/any_any_high"]
	fin_label = [r'$\alpha,\beta,\gamma = 0,300,0$'+' nm',r'$\alpha,\beta,\gamma = 30,300,30$'+' nm',r'$\alpha,\beta,\gamma = 30,300,10$'+' nm']
	marker = ["s","p","o","v","X"]
	hloc = ['left','left','right','right','right','right','right']
	zorder = [5,4,3,2,1,0]
	cf = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#CFECF9', '#7F7F7F', '#BCBD22', '#17BECF']
	norm = [3,3,3,3,2,2,2,2]
	yloc_r_ax1 = -2.6
	yloc_r_ax2 = -2.6
	xmin_ax1 = 1
	xmax_ax1 = 10**7 #5*10**4
	xmin_ax2 = 1
	xmax_ax2 = 10**7 #5*10**4
	xleg_ax1 = 3*10**3
	xleg_ax2 = 3*10**3
	ymin_ax2 = -14
	ymin_ax1 = -14
	ymax_ax1 = -2
	ymax_ax2 = -2
	flag = 1
elif (opt == '2'):
	itr = 'surf_high_v7'
	nxwindow = 10
	nxoverlap = 10
	fin_name = ["0_300/any_any_high","30_300_v2/30_50_10_300/any_any_high","30_300_v2/30_40_10_300/any_any_high","30_300_v2/30_10_300/any_any_high","30_300_v2/30_20_10_300/any_any_high"]
	fin_label = [r'$\alpha,\beta,\gamma,delta = 0,300,0,0$'+' nm',r'$\alpha,\beta,\gamma,delta = 30,300,10,50$'+' nm',r'$\alpha,\beta,\gamma,\delta = 30,300,10,40$'+' nm',r'$\alpha,\beta,\gamma,delta = 30,300,10,30$'+' nm',r'$\alpha,\beta,\gamma,\delta = 30,300,10,20$'+' nm']
	marker = ["s","p","o","v","X"]
	hloc = ['left','left','right','right','right',"right"]
	zorder = [5,4,3,2,1,0]
	cf = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#CFECF9', '#7F7F7F', '#BCBD22', '#17BECF']
	norm = [3,3,3,3,2,2,2,2]
	yloc_r_ax1 = -2.6
	yloc_r_ax2 = -2.6
	xmin_ax1 = 1
	xmax_ax1 = 10**7 #5*10**4
	xmin_ax2 = 1
	xmax_ax2 = 10**7 #5*10**4
	xleg_ax1 = 3*10**3
	xleg_ax2 = 3*10**3
	ymin_ax2 = -14
	ymin_ax1 = -14
	ymax_ax1 = -2
	ymax_ax2 = -2
	flag = 1
elif (opt == '3'):
	itr = 'surf_high_400'
	nxwindow = 10
	nxoverlap = 10
	fin_name = ["30_300_v2/0_0_0_400/any_any_high","30_300_v2/40_40_40_400/any_any_high","30_300_v2/40_40_10_400/any_any_high"]
	fin_label = [r'$\alpha,\beta,\gamma = 0,400,0$'+' nm',r'$\alpha,\beta,\gamma = 40,400,40$'+' nm',r'$\alpha,\beta,\gamma = 40,400,10$'+' nm']#,r'$\alpha,\beta,\gamma = 30,300,10$'+' nm']
	marker = ["s","p","o","v","X"]
	hloc = ['left','left','right','right','right',"right"]
	zorder = [5,4,3,2,1,0]
	cf = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#CFECF9', '#7F7F7F', '#BCBD22', '#17BECF']
	norm = [3,3,3,3,2,2,2,2]
	yloc_r_ax1 = -2.6
	yloc_r_ax2 = -2.6
	xmin_ax1 = 1
	xmax_ax1 = 10**7 #5*10**4
	xmin_ax2 = 1
	xmax_ax2 = 10**7 #5*10**4
	xleg_ax1 = 3*10**3
	xleg_ax2 = 3*10**3
	ymin_ax2 = -14
	ymin_ax1 = -14
	ymax_ax1 = -2
	ymax_ax2 = -2
	flag = 1
elif (opt == '4'):
	itr = 'surf_high_200'
	nxwindow = 10
	nxoverlap = 10
	fin_name = ["30_300_v2/0_0_0_200/any_any_high","30_300_v2/20_20_20_200/any_any_high","30_300_v2/20_20_5_200/any_any_high"]
	fin_label = [r'$\alpha,\beta,\gamma = 0,200,0$'+' nm',r'$\alpha,\beta,\gamma = 20,200,20$'+' nm',r'$\alpha,\beta,\gamma = 20,200,5$'+' nm']#,r'$\alpha,\beta,\gamma = 30,300,10$'+' nm']
	marker = ["s","p","o","v","X"]
	hloc = ['left','left','right','right','right',"right"]
	zorder = [5,4,3,2,1,0]
	cf = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#CFECF9', '#7F7F7F', '#BCBD22', '#17BECF']
	norm = [3,3,3,3,2,2,2,2]
	yloc_r_ax1 = -2.6
	yloc_r_ax2 = -2.6
	xmin_ax1 = 1
	xmax_ax1 = 10**7 #5*10**4
	xmin_ax2 = 1
	xmax_ax2 = 10**7 #5*10**4
	xleg_ax1 = 3*10**3
	xleg_ax2 = 3*10**3
	ymin_ax2 = -14
	ymin_ax1 = -14
	ymax_ax1 = -2
	ymax_ax2 = -2
	flag = 1
elif (opt == '5'):
	itr = 'surf_high_100'
	nxwindow = 10
	nxoverlap = 10
	fin_name = ["30_300_v2/0_0_0_100/any_any_high","30_300_v2/10_10_10_100/any_any_high","30_300_v2/10_10_3_100/any_any_high"]
	fin_label = [r'$\alpha,\beta,\gamma = 0,100,0$'+' nm',r'$\alpha,\beta,\gamma = 10,100,10$'+' nm',r'$\alpha,\beta,\gamma = 10,100,3$'+' nm']#,r'$\alpha,\beta,\gamma = 30,300,10$'+' nm']
	marker = ["s","p","o","v","X"]
	hloc = ['left','left','right','right','right',"right"]
	zorder = [5,4,3,2,1,0]
	cf = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#CFECF9', '#7F7F7F', '#BCBD22', '#17BECF']
	norm = [3,3,3,3,2,2,2,2]
	yloc_r_ax1 = -2.6
	yloc_r_ax2 = -2.6
	xmin_ax1 = 1
	xmax_ax1 = 10**7 #5*10**4
	xmin_ax2 = 1
	xmax_ax2 = 10**7 #5*10**4
	xleg_ax1 = 3*10**3
	xleg_ax2 = 3*10**3
	ymin_ax2 = -14
	ymin_ax1 = -14
	ymax_ax1 = -2
	ymax_ax2 = -2
	flag = 1
elif (opt == '6'):
	itr = 'surf_high_150'
	nxwindow = 10
	nxoverlap = 10
	fin_name = ["30_300_v2/0_0_0_150/any_any_high","30_300_v2/15_15_15_150/any_any_high","30_300_v2/15_15_5_150/any_any_high"]
	fin_label = [r'$\alpha,\beta,\gamma = 0,150,0$'+' nm',r'$\alpha,\beta,\gamma = 15,150,15$'+' nm',r'$\alpha,\beta,\gamma = 15,150,5$'+' nm']#,r'$\alpha,\beta,\gamma = 30,300,10$'+' nm']
	marker = ["s","p","o","v","X"]
	hloc = ['left','left','right','right','right',"right"]
	zorder = [5,4,3,2,1,0]
	cf = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#CFECF9', '#7F7F7F', '#BCBD22', '#17BECF']
	norm = [3,3,3,3,2,2,2,2]
	yloc_r_ax1 = -2.6
	yloc_r_ax2 = -2.6
	xmin_ax1 = 1
	xmax_ax1 = 10**7 #5*10**4
	xmin_ax2 = 1
	xmax_ax2 = 10**7 #5*10**4
	xleg_ax1 = 3*10**3
	xleg_ax2 = 3*10**3
	ymin_ax2 = -14
	ymin_ax1 = -14
	ymax_ax1 = -2
	ymax_ax2 = -2
	flag = 1
elif (opt == '7'):
	itr = 'surf_high_250'
	nxwindow = 10
	nxoverlap = 10
	fin_name = ["30_300_v2/0_0_0_250/any_any_high","30_300_v2/25_25_7_250/any_any_high","30_300_v2/25_25_8_250/any_any_high"]
	fin_label = [r'$\alpha,\beta,\gamma = 0,250,0$'+' nm',r'$\alpha,\beta,\gamma = 25,250,7$'+' nm',r'$\alpha,\beta,\gamma = 25,250,8$'+' nm']#,r'$\alpha,\beta,\gamma = 30,300,10$'+' nm']
	marker = ["s","p","o","v","X"]
	hloc = ['left','left','right','right','right',"right"]
	zorder = [5,4,3,2,1,0]
	cf = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#CFECF9', '#7F7F7F', '#BCBD22', '#17BECF']
	norm = [3,3,3,3,2,2,2,2]
	yloc_r_ax1 = -2.6
	yloc_r_ax2 = -2.6
	xmin_ax1 = 1
	xmax_ax1 = 10**7 #5*10**4
	xmin_ax2 = 1
	xmax_ax2 = 10**7 #5*10**4
	xleg_ax1 = 3*10**3
	xleg_ax2 = 3*10**3
	ymin_ax2 = -14
	ymin_ax1 = -14
	ymax_ax1 = -2
	ymax_ax2 = -2
	flag = 1
elif (opt == '500'):
	itr = 'surf_high_500'
	nxwindow = 10
	nxoverlap = 10
	fin_name = ["30_300_v2/0_0_0_500/any_any_high","30_300_v2/50_50_50_500/any_any_high","30_300_v2/50_50_25_500/any_any_high","30_300_v2/50_50_16_500/any_any_high","30_300_v2/50_50_10_500/any_any_high","30_300_v2/50_50_6_500/any_any_high"]#,"30_300_v2/25_25_8_250/any_any_high"]
	fin_label = [r'$\alpha,\beta,\gamma = 0,500,0$'+' nm',r'$\alpha,\beta,\gamma = 50,500,50$'+' nm',r'$\alpha,\beta,\gamma = 50,500,25$'+' nm',r'$\alpha,\beta,\gamma = 50,500,16$'+' nm',r'$\alpha,\beta,\gamma = 50,500,10$'+' nm',r'$\alpha,\beta,\gamma = 50,500,6$'+' nm']#,r'$\alpha,\beta,\gamma = 25,250,8$'+' nm']#,r'$\alpha,\beta,\gamma = 30,300,10$'+' nm']
	marker = ["s","p","o","v","X","."]
	hloc = ['left','left','right','right','right',"right","right"]
	zorder = [6,5,4,3,2,1,0]
	cf = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#CFECF9', '#7F7F7F', '#BCBD22', '#17BECF']
	norm = [3,3,3,3,2,2,2,2]
	yloc_r_ax1 = -2.6
	yloc_r_ax2 = -2.6
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
elif (opt == '500_0_10'):
	itr = 'surf_high_'+opt
	nxwindow = 10
	nxoverlap = 10
	fin_name = ["30_300_v2/0_0_0_500/any_any_high","30_300_v2/50_50_10_500/any_any_high"]#,"30_300_v2/25_25_8_250/any_any_high"]
	#fin_label = [r'$\alpha,\beta,\gamma = 0,500,0$'+' nm',r'$\alpha,\beta,\gamma = 50,500,10$'+' nm']#,r'$\alpha,\beta,\gamma = 25,250,8$'+' nm']#,r'$\alpha,\beta,\gamma = 30,300,10$'+' nm']
	fin_label = [r'$Smooth$',r'$Finned$']
	marker = ["s","p","o","v","X","."]
	hloc = ['left','right','right','right','right',"right","right"]
	zorder = [6,5,4,3,2,1,0]
	cf = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#CFECF9', '#7F7F7F', '#BCBD22', '#17BECF']
	norm = [3,3,3,3,2,2,2,2]
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
elif (opt == '500_0_4'):
	itr = 'surf_high_'+opt
	nxwindow = 10
	nxoverlap = 10
	fin_name = ["30_300_v2/0_0_0_500/any_any_high","30_300_v2/50_50_4_500/any_any_high"]#,"30_300_v2/25_25_8_250/any_any_high"]
	fin_label = [r'$\alpha,\beta,\gamma = 0,500,0$'+' nm',r'$\alpha,\beta,\gamma = 50,500,4$'+' nm']#,r'$\alpha,\beta,\gamma = 25,250,8$'+' nm']#,r'$\alpha,\beta,\gamma = 30,300,10$'+' nm']
	marker = ["s","p","o","v","X","."]
	hloc = ['left','right','right','right','right',"right","right"]
	zorder = [6,5,4,3,2,1,0]
	cf = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#CFECF9', '#7F7F7F', '#BCBD22', '#17BECF']
	norm = [3,3,3,3,2,2,2,2]
	yloc_r_ax1 = -2.6
	yloc_r_ax2 = -2.6
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
elif (opt == '500_0_4_10_25_50'):
	itr = 'surf_high_'+opt
	nxwindow = 10
	nxoverlap = 10
	fin_name = ["30_300_v2/0_0_0_500/any_any_high","30_300_v2/50_50_4_500/any_any_high","30_300_v2/50_50_10_500/any_any_high","30_300_v2/50_50_25_500/any_any_high","30_300_v2/50_50_50_500/any_any_high"]
	fin_label = [r'$\alpha,\beta,\gamma = 0,500,0$'+' nm',r'$\alpha,\beta,\gamma = 50,500,4$'+' nm',r'$\alpha,\beta,\gamma = 50,500,10$'+' nm',r'$\alpha,\beta,\gamma = 50,500,25$'+' nm',r'$\alpha,\beta,\gamma = 50,500,50$'+' nm']
	marker = ["s","p","o","v","X","."]
	hloc = ['left','right','right','right','right',"right","right"]
	zorder = [6,5,4,3,2,1,0]
	cf = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#CFECF9', '#7F7F7F', '#BCBD22', '#17BECF']
	norm = [3,3,3,3,2,2,2,2]
	yloc_r_ax1 = -2.6
	yloc_r_ax2 = -2.6
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
elif (opt == '500_0new_4_10_25_50'):
	itr = 'surf_high_'+opt
	fin_name = ["30_300_v2/0_0_0_500new/any_any_high","30_300_v2/50_50_4_500/any_any_high","30_300_v2/50_50_10_500/any_any_high","30_300_v2/50_50_25_500/any_any_high","30_300_v2/50_50_50_500/any_any_high"]
	fin_label = [r'$\alpha,\beta,\gamma = 0,500,0$'+' nm',r'$\alpha,\beta,\gamma = 50,500,4$'+' nm',r'$\alpha,\beta,\gamma = 50,500,10$'+' nm',r'$\alpha,\beta,\gamma = 50,500,25$'+' nm',r'$\alpha,\beta,\gamma = 50,500,50$'+' nm']
	flag = 1


# In[3]:


opt = sys.argv[1]
if (opt == '500_0new_4_10_25_50'):
	itr = 'surf_high_'+opt
	fin_name = ["30_300_v2/0_0_0_500new/any_any_high","30_300_v2/50_50_4_500/any_any_high","30_300_v2/50_50_10_500/any_any_high","30_300_v2/50_50_25_500/any_any_high","30_300_v2/50_50_50_500/any_any_high"]
	fin_label = [r'$\alpha,\beta,\gamma = 0,500,0$'+' nm',r'$\alpha,\beta,\gamma = 50,500,4$'+' nm',r'$\alpha,\beta,\gamma = 50,500,10$'+' nm',r'$\alpha,\beta,\gamma = 50,500,25$'+' nm',r'$\alpha,\beta,\gamma = 50,500,50$'+' nm']
	flag = 1


# In[ ]:


def getfiles(dirname):
	names = glob.glob(dirname+"/traj*_*")
	return names


# In[ ]:


def readfile(name,flag):
	f = open(name,'r')
	i = 0
	flag = 0
	pathtime = np.nan
	lasttime = np.nan
	firsttime = np.nan
	lines = reversed(f.readlines())
	for line in lines:
		i = i+1
		if (i==1 and line=="Finished"):
			flag = 1
		secondlastline = "0"
		if (i==3 and flag==1):
			secondlastline = line
			if (len(secondlastline.split("\t"))==5):
				lasttime = secondlastline.split("\t")[0]
	f.close()
	f = open(name,'r')
	if flag == 1:
		lines = f.readlines()
		i = 0
		for line in lines:
			i = i+1
			if(i==2 and flag==1):
				secondline = line
				firsttime = secondline.split("\t")[0]
	pathtime = float(lasttime)-float(firsttime)
	return float(pathtime)


# In[ ]:


def escape_time_cumuprob(dirname,flag):
	names = getfiles(dirname)
	t = []
	for name in names:
		t+=[readfile(name,flag)]
	t = np.array(t)
	t = t[np.isfinite(t)]
	return t


# In[ ]:


def overlaphistogram(xy,xwindow,xoverlap,dx):
	xedge = np.exp(np.arange(np.log(0.1),np.log(xy.max())+1.2,0.4))
	ynew,xedge = np.histogram(xy,bins=xedge)
	xcenter = (xedge[1:]+xedge[:-1])*0.5
	del xedge
	return xcenter,ynew


# In[ ]:


def plotcumtime(dirname,lab):
	global itr,fin_name,fin_label,nofin_name,nofin_label,flag
	xy = escape_time_cumuprob(dirname,flag)
	xy = xy[:min(144000,len(xy))]
	avgt = np.average(xy)
	avgterr = np.sqrt(np.var(xy)/len(xy))
	print dirname.split('/')[-2],"{:.3E}".format(avgt),'+-',"{:.3E}".format(avgterr)
	print "For ",len(xy),' : ',min(xy),max(xy)


# In[ ]:


for i in range(len(fin_name)):
	plotcumtime(dirname=fin_name[i],lab=fin_label[i])

