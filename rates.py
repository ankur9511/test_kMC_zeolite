#!/usr/bin/env python2
# coding: utf-8

# In[ ]:


import numpy as np
import sys


# In[ ]:


flag = sys.argv[1]


# In[ ]:


def Forester_Smith_High_des():
	ratelist = {}
	ratelist['ZI'] = {}
	ratelist['ZZ'] = {}
	ratelist['ZS'] = {}
	ratelist['ZGI']= {}
	ratelist['GI'] = {}
	ratelist['ZI']['ZZ'] = 1.1068*10**6
	ratelist['ZZ']['ZI'] = 1.5834*10**8
	ratelist['ZI']['ZS'] = 2.3067*10**6
	ratelist['ZS']['ZI']  = 2.7333*10**8
	##### Ads Des rate
	ratelist['ZI']['ZGI'] = 10**9
	ratelist['ZGI']['ZI'] = 0
	ratelist['ZGI']['ZGI'] = {}
	ratelist['ZGI']['ZGI']['xy'] = 0
	ratelist['ZGI']['ZGI']['z']  = 0
	ratelist['ZGI']['GI'] = {}
	ratelist['ZGI']['GI']['xy']  = 0
	ratelist['ZGI']['GI']['z']   = 0
	ratelist['GI']['ZGI'] = {}
	ratelist['GI']['ZGI']['xy']  = 0
	ratelist['GI']['ZGI']['z']   = 0
	ratelist['GI']['GI']         = {}
	ratelist['GI']['GI']['xy']   = 0
	ratelist['GI']['GI']['z']    = 0
	np.save("FSHdes_ratelist",ratelist)


# In[ ]:


def Forester_Smith_Low_des():
	ratelist = {}
	ratelist['ZI'] = {}
	ratelist['ZZ'] = {}
	ratelist['ZS'] = {}
	ratelist['ZGI']= {}
	ratelist['GI'] = {}
	ratelist['ZI']['ZZ'] = 1.1068*10**6
	ratelist['ZZ']['ZI'] = 1.5834*10**8
	ratelist['ZI']['ZS'] = 2.3067*10**6
	ratelist['ZS']['ZI']  = 2.7333*10**8
	##### Ads Des rate
	ratelist['ZI']['ZGI'] = 10**5
	ratelist['ZGI']['ZI'] = 0
	ratelist['ZGI']['ZGI'] = {}
	ratelist['ZGI']['ZGI']['xy'] = 0
	ratelist['ZGI']['ZGI']['z']  = 0
	ratelist['ZGI']['GI'] = {}
	ratelist['ZGI']['GI']['xy']  = 0
	ratelist['ZGI']['GI']['z']   = 0
	ratelist['GI']['ZGI'] = {}
	ratelist['GI']['ZGI']['xy']  = 0
	ratelist['GI']['ZGI']['z']   = 0
	ratelist['GI']['GI']         = {}
	ratelist['GI']['GI']['xy']   = 0
	ratelist['GI']['GI']['z']    = 0
	np.save("FSLdes_ratelist",ratelist)
	
if (flag == "HighDesorption"):
	Forester_Smith_High_des()
elif (flag == "LowDesorption"):
	Forester_Smith_Low_des()

