#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 21:06:50 2017

@author: eweng
"""
import numpy as np
import csv
import matplotlib.pyplot as plt

f=open('../output/Ecosystem_yearly.csv')
next(f) # skip headings
yearlydata=list(csv.reader(f,delimiter=','))
f.close()

#with open('../output/Ecosystem_yearly.csv') as csvfile:
#    readCSV = csv.reader(csvfile, delimiter=',')
## data transformation

ecodatayr=np.array(yearlydata)
xyear = ecodatayr[:,0].astype(np.float)
GPP = ecodatayr[:,3].astype(np.float)
NPP = ecodatayr[:,4].astype(np.float)
Rh  = ecodatayr[:,5].astype(np.float)
Transp = ecodatayr[:,8].astype(np.float)
Evap   = ecodatayr[:,9].astype(np.float)
Runoff = ecodatayr[:,10].astype(np.float)
#%% Plot

font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }
#x = np.linspace(0.0, 5.0, 100)
#y = np.cos(2 * np.pi * x) * np.exp(-x)
plt.figure(1)
plt.plot(xyear, GPP, 'b',xyear, NPP, 'r',xyear, Rh, 'k')
plt.legend( ('GPP', 'NPP', 'Rh') )
plt.title('C fluxes', fontdict=font)
plt.xlabel('Year', fontdict=font)
plt.ylabel('C flux (kgC $m^{-2} yr^{-1}$)', fontdict=font)
# Tweak spacing to prevent clipping of ylabel
#plt.subplots_adjust(left=0.15)
plt.show()

plt.figure(2)
plt.plot(xyear, Transp, 'b*', xyear, Evap, 'ro',xyear, Runoff, 'k+')
plt.legend( ('Transp', 'Evap', 'Runoff') )
plt.title('Water fluxes', fontdict=font)
plt.xlabel('Year', fontdict=font)
plt.ylabel('Water flux (mm $yr^{-1}$)', fontdict=font)
plt.show()
