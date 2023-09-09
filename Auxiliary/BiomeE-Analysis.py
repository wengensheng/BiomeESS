#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 06/17/2023
  for BiomeE output
@author: eweng
"""
import numpy as np
import csv
import matplotlib.pyplot as plt
from cycler import cycler

linestyle_cycler = cycler('linestyle',['-','-.','--',':'])

font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 14,
        }

#%% Path and file names
fpath = '../output/'
#fpath = './BiomeESimulations/'
fout  = '../output/'
dtype = ['Ecosystem_yearly','Cohort_yearly']
N_pfts = 4 # total PFTs at one site, 4
N_Layers = 3
PI = 3.1415926
DBHbins=[0.0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.5,2.0,99]

# Site and experiments
siteID = 'BCI_hydro' #'ORNL_test'
PFTID = ['PFT1','PFT2','PFT3','PFT4']

#%% Retrieve data
# File name
expID = siteID
fname1 = fpath + expID + '_' + dtype[0] + '.csv'
fname2 = fpath + expID + '_' + dtype[1] + '.csv'


# Land data yearly
f=open(fname1)
next(f) # skip headings
LandYrV=(np.array(list(csv.reader(f,delimiter=','))))
f.close()

# Total cohort yearly
f=open(fname2)
next(f) # skip headings
CCYrV=(np.array(list(csv.reader(f,delimiter=','))))
f.close()

# Land Data dimensions
rows = len(LandYrV)
col  = len(LandYrV[1]) - 1
LandYr = LandYrV[0:rows,1:col].astype(np.float)
totYrs = rows

# Total cohort yearly dimensions
rows = len(CCYrV)
col  = len(CCYrV[1]) - 1
CCYr = CCYrV[0:rows,0:col].astype(np.float)
totCCL = rows

#% Data calculation and analysis
LandC = np.copy(LandYr[:,[5,6,7,14,15]]) # GPP, Rauto, Rh, PlantC, SoilC

AGB = np.zeros(totYrs)
Rh  = np.zeros(totYrs)
nbp = np.zeros(totYrs)

BA = np.zeros((totYrs, N_pfts))
CA = np.zeros((totYrs, N_pfts))
LAI = np.zeros((totYrs, N_pfts))
height = np.zeros((totYrs, N_pfts))
cveg = np.zeros((totYrs, N_pfts))
cwood = np.zeros((totYrs, N_pfts))
gpp = np.zeros((totYrs, N_pfts))
npp = np.zeros((totYrs, N_pfts))
WBgrowth = np.zeros((totYrs, N_pfts))
BAgrowth = np.zeros((totYrs, N_pfts))
WDmort = np.zeros((totYrs, N_pfts))
DENmort = np.zeros((totYrs, N_pfts))

Nstem1 = np.zeros((totYrs, N_pfts))
Nstem2 = np.zeros((totYrs, N_pfts))

CA_grow = np.zeros((totYrs, N_Layers))
CA_mort = np.zeros((totYrs, N_Layers))
DEN_L  = np.zeros((totYrs, N_Layers))

cwood_size = np.zeros((totYrs,16))
nstem_size = np.zeros((totYrs,16))
cmort_size = np.zeros((totYrs,16))

stemmort_size = np.zeros((totYrs,16, N_pfts))

for i in range(totCCL):
    iYr  = int(CCYr[i,1])-1
    iPFT = int(CCYr[i,4])-1
    iLayer = min(2,int(CCYr[i,5])-1)
    Nstem1[iYr,iPFT] = Nstem1[iYr,iPFT]  + CCYr[i,6]
    BA[iYr,iPFT]  = BA[iYr,iPFT]  + CCYr[i,6]*PI*0.25*CCYr[i,11]**2
    height[iYr,iPFT] = max(height[iYr,iPFT],CCYr[i,12])
    CA[iYr,iPFT]  = CA[iYr,iPFT]  + CCYr[i,6]*CCYr[i,13]/10000
    LAI[iYr,iPFT] = LAI[iYr,iPFT] + CCYr[i,6]*CCYr[i,14]/10000
    gpp[iYr,iPFT] = gpp[iYr,iPFT] + CCYr[i,6]*CCYr[i,22]/10000
    npp[iYr,iPFT] = npp[iYr,iPFT] + CCYr[i,6]*CCYr[i,23]/10000
    cwood[iYr,iPFT] = cwood[iYr,iPFT] + CCYr[i,6]*(CCYr[i,17]+CCYr[i,18])/10000
    cveg[iYr,iPFT]  = cveg[iYr,iPFT]  + CCYr[i,6]*np.sum(CCYr[i,15:21])/10000
    WBgrowth[iYr,iPFT] = WBgrowth[iYr,iPFT] + CCYr[i,6]* CCYr[i,24] * CCYr[i,28]/10000
    BAgrowth[iYr,iPFT] = BAgrowth[iYr,iPFT] + CCYr[i,6]* CCYr[i,9] #/10000.0
    WDmort[iYr,iPFT] = WDmort[iYr,iPFT] + CCYr[i,6]*(CCYr[i,17]+CCYr[i,18])*CCYr[i,29]/10000.
    DENmort[iYr,iPFT] = DENmort[iYr,iPFT] + CCYr[i,6] * CCYr[i,29]/10000.
    CA_grow[iYr,iLayer] = CA_grow[iYr,iLayer] + CCYr[i,6] * (1.0-CCYr[i,29])*CCYr[i,10]
    CA_mort[iYr,iLayer] = CA_mort[iYr,iLayer] + CCYr[i,6] * CCYr[i,29]*CCYr[i,13]
    DEN_L[iYr,iLayer] = DEN_L[iYr,iLayer]  + CCYr[i,6]
    for j in range(len(DBHbins)-1):
        if CCYr[i,11]>=DBHbins[j] and CCYr[i,11]<DBHbins[j+1]:
            cwood_size[iYr,j] = cwood_size[iYr,j] + CCYr[i,6]*(CCYr[i,17]+CCYr[i,18])/10000.0 # KgC/m2
            nstem_size[iYr,j] = nstem_size[iYr,j] + CCYr[i,6]
            cmort_size[iYr,j] = cmort_size[iYr,j] + CCYr[i,6]*(CCYr[i,17]+CCYr[i,18])*CCYr[i,29]/10000. # KgC/m2/yr
            stemmort_size[iYr,j,iPFT] = stemmort_size[iYr,j,iPFT] + CCYr[i,5] * CCYr[i,29]/10000. # stems/yr
            break

AGB[:] = np.sum(cveg[:,1:4],axis=1)*0.75
Rh[:]  = LandC[:,2]
nbp[:] = LandC[:,0]-LandC[:,1]-LandC[:,2] # NBP = GPP - NPP - Rh

meanWDCclass = np.mean(cwood_size[totYrs-60:totYrs,:],axis=0)
meanDen = np.mean(nstem_size[totYrs-60:totYrs,:],axis=0)
meanDen[0:2] = 0.0

#%% Plot
xyear = np.arange(totYrs)

# Vegetation structure
plt.figure(1) #
plt.clf()
plt.subplot(221)
plt.plot(xyear, BA)
plt.xlabel('Year', fontdict=font)
plt.ylabel('Basal area\n (m$^{2}$ ha$^{-1}$)', fontdict=font)

plt.subplot(222)
plt.plot(xyear, CA)
plt.legend((PFTID),loc=0,ncol=1)
plt.xlabel('Year', fontdict=font)
plt.ylabel('Crown area\n (m$^{2}$ m$^{-2}$)', fontdict=font)

plt.subplot(223)
plt.plot(xyear, LAI)
#plt.legend((PFTID),loc=0,ncol=1)
plt.xlabel('Year', fontdict=font)
plt.ylabel('Leaf area\n (m$^{2}$ m$^{-2}$)', fontdict=font)

plt.subplot(224)
plt.plot(xyear, height)
#plt.legend((PFTID),loc=0,ncol=1)
plt.xlabel('Year', fontdict=font)
plt.ylabel('Height (95%) (m)', fontdict=font)

# Vegetation carbon
plt.figure(2) #
plt.clf()
plt.subplot(221)
plt.plot(xyear,gpp)
plt.xlabel('Year', fontdict=font)
plt.ylabel('GPP\n (KgC m$^{-2}$ yr$^{-1}$)', fontdict=font)

plt.subplot(222)
plt.plot(xyear, npp)
plt.legend((PFTID),loc=0,ncol=1)
plt.xlabel('Year', fontdict=font)
plt.ylabel('NPP\n (KgC m$^{-2}$ yr$^{-1}$)', fontdict=font)

plt.subplot(223)
plt.plot(xyear, WBgrowth)
#plt.legend((PFTID),loc=0,ncol=1)
plt.xlabel('Year', fontdict=font)
plt.ylabel('Wood growth \n(KgC m$^{2}$ yr$^{-1}$)', fontdict=font)

plt.subplot(224)
plt.plot(xyear, cveg)
#plt.legend((PFTID),loc=0,ncol=1)
plt.xlabel('Year', fontdict=font)
plt.ylabel('Biomass\n (KgC m$^{-2}$)', fontdict=font)

# Ecosystem C
plt.figure(3) #
plt.clf()
plt.subplot(221)
plt.plot(xyear, LandC[:,0])
plt.xlabel('Year', fontdict=font)
plt.ylabel('GPP\n(KgC m$^{-2}$ yr$^{-1}$)', fontdict=font)

plt.subplot(222)
plt.plot(xyear, LandC[:,1])
plt.xlabel('Year', fontdict=font)
plt.ylabel('Rauto\n(KgC m$^{-2}$ yr$^{-1}$)', fontdict=font)

plt.subplot(223)
plt.plot(xyear, LandC[:,3])
plt.xlabel('Year', fontdict=font)
plt.ylabel('Plant C \n(KgC m$^{-2}$ yr$^{-1}$)', fontdict=font)

plt.subplot(224)
plt.plot(xyear, LandC[:,4])
plt.xlabel('Year', fontdict=font)
plt.ylabel('Soil C \n(KgC m$^{-2}$ yr$^{-1}$)', fontdict=font)


# Size classes
xbins = np.arange(16)
plt.figure(4) #
plt.clf()
plt.subplot(221)
plt.scatter(xbins, meanDen)
#plt.yscale("log")
plt.xlabel('Bin', fontdict=font)
plt.ylabel('Stems/ha', fontdict=font)

plt.subplot(222)
plt.scatter(xbins, meanWDCclass)
plt.xlabel('Bin', fontdict=font)
plt.ylabel('Wood C (KgC m$^{-2}$)', fontdict=font)

plt.subplot(223)
plt.plot(xyear, CA_grow[:,0], CA_mort[:,0])
plt.xlabel('Year', fontdict=font)
plt.ylabel('dCA (m$^{2}$/ha)', fontdict=font)
plt.ylim((0,500))

plt.subplot(224)
plt.plot(xyear, CA_grow[:,1], CA_mort[:,1])
plt.xlabel('Year', fontdict=font)
plt.ylabel('dCA (m$^{2}$/ha)', fontdict=font)

# Growth and mortality
plt.figure(5) #
plt.clf()
plt.subplot(221)
plt.plot(xyear, WBgrowth)
#plt.legend((PFTID),loc=0,ncol=1)
plt.xlabel('Year', fontdict=font)
plt.ylabel('Wood growth \n(KgC m$^{2}$ yr$^{-1}$)', fontdict=font)

BAgrowth[:,0] = 0.0
plt.subplot(222)
plt.plot(xyear,BAgrowth)
plt.xlabel('Year', fontdict=font)
plt.ylabel('dBA\n (m$^{2}$ ha$^{-1}$ yr$^{-1}$)', fontdict=font)

plt.subplot(223)
plt.plot(xyear, WDmort)
plt.legend((PFTID),loc=0,ncol=1)
plt.xlabel('Year', fontdict=font)
plt.ylabel('WD mort\n (KgC m$^{-2}$ yr$^{-1}$)', fontdict=font)

DENmort[:,0] = 0.0
plt.subplot(224)
plt.plot(xyear, DENmort)
#plt.legend((PFTID),loc=0,ncol=1)
plt.xlabel('Year', fontdict=font)
plt.ylabel('Den mort\n (n m$^{-2}$ yr$^{-1}$)', fontdict=font)

# Size class
plt.figure(6) #
plt.clf()

plt.subplot(121)
plt.plot(xyear, nstem_size)
plt.yscale("log")
plt.xlabel('Year', fontdict=font)
plt.ylabel('Stems/ha', fontdict=font)

plt.subplot(122)
plt.plot(xyear, cwood_size)
plt.xlabel('Year', fontdict=font)
plt.ylabel('Wood C (KgC m$^{-2}$)', fontdict=font)


# Crown dynamics
plt.figure(7) #
plt.clf()
plt.subplot(221)
plt.plot(xyear, CA[:,2])
plt.xlabel('Year', fontdict=font)
plt.ylabel('Crown area (m$^{2}$/ha)', fontdict=font)

plt.subplot(222)
plt.plot(xyear, DEN_L[:,0:2])
plt.legend(('1st layer','2nd layer'),loc=0,ncol=1)
plt.yscale('log')
plt.xlabel('Year', fontdict=font)
plt.ylabel('Density (Stems/ha)', fontdict=font)

plt.subplot(223)
plt.plot(xyear, CA_grow[:,0], CA_mort[:,0])
plt.legend(('CA_grow','CA_mort'),loc=0,ncol=1)
#plt.yscale('log')
plt.xlabel('Year', fontdict=font)
plt.ylabel('dCA-Layer1 (m$^{2}$/ha)', fontdict=font)
plt.ylim((0,400))

plt.subplot(224)
plt.plot(xyear, CA_grow[:,1], CA_mort[:,1])
plt.legend(('CA_grow','CA_mort'),loc=0,ncol=1)
plt.xlabel('Year', fontdict=font)
plt.ylabel('dCA-Layer2 (m$^{2}$/ha)', fontdict=font)

