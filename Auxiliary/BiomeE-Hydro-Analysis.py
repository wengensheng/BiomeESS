#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 07/21/2022
  for BiomeE plant hydraulics
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

ls_cycler = cycler('linestyle',
            [(0,()), # solid
            (0, (1, 10)), # loosely dotted
            (0, (1, 5)), # dotted
            (0, (1, 1)), # densely dotted
            (0, (5, 10)), # loosely dashed
            (0, (5, 5)), # dashed
            (0, (5, 1)), # densely dashed
            (0, (3, 10, 1, 10)), # loosely dashdotted
            (0, (3, 5, 1, 5)), # dashdotted
            (0, (3, 1, 1, 1)), # densely dashdotted
            (0, (3, 10, 1, 10, 1, 10)), # loosely dashdotdotted
            (0, (3, 5, 1, 5, 1, 5)), # dashdotdotted
            (0, (3, 1, 1, 1, 1, 1))] # densely dashdotdotted 
                  )
#%% Path and file names
fpath = '../output/'
expID = ['PFT3-iWTC0iKx0','PFT3-cWTC0iKx0','5PFTs-iWTC0iKx0','Single5PFTs-iW0iK0-2']
fname = ['Summary','Patch_yearly','FirstCohortYearly','FirstCohortDaily', \
         'FirstCohortHourly','Cohort_yearly_clean']
runID = '_BCI_hydro'

PFTID = ['WD200','WD250','WD300','WD350','WD400']

N_tests = 3
#yrmax = 612
#EcoData = np.zeros((N_tests,yrmax,42)) # 42 is "summary" columns
# Read in data

# Patch data yearly
fname1 = fpath+expID[0]+'/'+fname[1]+runID+'.csv'
f=open(fname1)
next(f) # skip headings
PatchYrV=(np.array(list(csv.reader(f,delimiter=','))))
f.close()

fname1 = fpath+expID[1]+'/'+fname[1]+runID+'.csv'
f=open(fname1)
next(f) # skip headings
PatchYrC=(np.array(list(csv.reader(f,delimiter=','))))
f.close()

# First cohort yearly
fname1 = fpath+expID[0]+'/'+fname[2]+runID+'.csv'
f=open(fname1)
next(f) # skip headings
CCYrV=(np.array(list(csv.reader(f,delimiter=','))))
f.close()

fname1 = fpath+expID[1]+'/'+fname[2]+runID+'.csv'
f=open(fname1)
next(f) # skip headings
CCYrC=(np.array(list(csv.reader(f,delimiter=','))))
f.close()

# First cohort hourly
fname1 = fpath+expID[0]+'/'+fname[4]+runID+'.csv'
f=open(fname1)
next(f) # skip headings
CCHrV=(np.array(list(csv.reader(f,delimiter=','))))
f.close()

fname1 = fpath+expID[1]+'/'+fname[4]+runID+'.csv'
f=open(fname1)
next(f) # skip headings
CCHrC=(np.array(list(csv.reader(f,delimiter=','))))
f.close()

#%% Retrieve data

# Patch Data dimensions
rows = len(PatchYrV)
col  = len(PatchYrV[1]) - 1
PatchYrVwtc0=PatchYrV[0:rows,0:col].astype(np.float)
PatchYrCwtc0=PatchYrC[0:rows,0:col].astype(np.float)

# First cohort yearly dimensions
rows = len(CCYrV)
col  = len(CCYrV[1]) - 1
CCYrVwtc0 = CCYrV[0:rows,0:25].astype(np.float)
CCYrCwtc0 = CCYrC[0:rows,0:25].astype(np.float)


#%% Plotting succession at increasing mortality rate and constant mort.

fname1 = fpath+'5PFTs-iWTC0iKx0'+'/'+fname[0]+runID+'.csv'
f=open(fname1)
next(f) # skip headings
RawData=(np.array(list(csv.reader(f,delimiter=','))))
f.close()
# Data dimensions
rows = len(RawData)
col  = len(RawData[1]) - 1
yrmax = rows

# Get data to new arrays
xyear = RawData[0:yrmax,0].astype(np.float)
EcoData = RawData[0:yrmax,0:col].astype(np.float)

fname1 = fpath+'5PFTs-iW0iK0-SameU'+'/'+fname[0]+runID+'.csv'
f=open(fname1)
next(f) # skip headings
RawData=(np.array(list(csv.reader(f,delimiter=','))))
f.close()
# Data dimensions
rows = len(RawData)
col  = len(RawData[1]) - 1
yrmax = rows

# Get data to new arrays
SameUData = RawData[0:yrmax,0:col].astype(np.float)

#%% Plot successional patterns
plt.figure(1) #

plt.rc('axes', prop_cycle=(cycler('color', ['b', 'y', 'g', 'm', 'r'])+
    cycler('linestyle', ['-', '-.', '--', ':',(0, (3, 1, 1, 1, 1, 1))])))
plt.clf()


plt.subplot(211)
plt.plot(xyear, SameUData[:,1:6])
plt.legend((PFTID),title='PFT',loc=0,ncol=1) #
#plt.ylim((0.0,60))
plt.title('a', x=0.1, y=0.85,fontdict=font)
plt.tick_params(labelbottom='off')
#plt.xlabel('Year', fontdict=font)
plt.ylabel('Basal area\n(m$^{2}$ ha$^{-1}$)', fontdict=font)

plt.subplot(212)
plt.plot(xyear, EcoData[:,1:6])
#plt.legend((PFTID),title='PFT',loc=0,ncol=1) #
#plt.ylim((0.0,60))
plt.title('b', x=0.1, y=0.85,fontdict=font)
#plt.tick_params(labelbottom='off')
plt.xlabel('Year', fontdict=font)
plt.ylabel('Basal area\n(m$^{2}$ ha$^{-1}$)', fontdict=font)


#%% Plotting yearly patch
plt.figure(2) #
#plt.rc('axes', prop_cycle=(cycler('color', ['b','r', 'g',  'y', 'm'])+
#    cycler('linestyle', ['-', '--', ':', '-.',(0, (3, 1, 1, 1, 1, 1))])))
plt.clf()

plt.subplot(221)
plt.plot(xyear, PatchYrVwtc0[:,5],PatchYrCwtc0[:,5])
plt.legend((['Increasing WTC0', 'Constant WTC0']),loc=0,ncol=1) # title='PFT',
plt.ylim((0.0,4))
plt.title('a', x=0.1, y=0.85,fontdict=font)
plt.tick_params(labelbottom='off')
#plt.xlabel('Year', fontdict=font)
plt.ylabel('GPP\n(Kg C m$^{-2}$ yr$^{-1}$)', fontdict=font)

plt.subplot(222)
plt.plot(xyear, PatchYrVwtc0[:,11],PatchYrCwtc0[:,11])
#plt.legend((expID[0:2]),loc=0,ncol=1) # title='PFT',
plt.ylim((0.0,1250))
plt.title('b', x=0.1, y=0.85,fontdict=font)
plt.tick_params(labelbottom='off')
#plt.xlabel('Year', fontdict=font)
plt.ylabel('Transp (mm yr$^{-1}$)', fontdict=font)

plt.subplot(223)
plt.plot(xyear, PatchYrVwtc0[:,14],PatchYrCwtc0[:,14])
#plt.legend((expID[0:2]),loc=0,ncol=1) # title='PFT',
#plt.ylim((0.0,60))
plt.title('c', x=0.1, y=0.85,fontdict=font)
#plt.tick_params(labelbottom='off')
plt.xlabel('Year', fontdict=font)
plt.ylabel('Plant C (Kg C m$^{-2}$)', fontdict=font)

plt.subplot(224)
plt.plot(xyear, PatchYrVwtc0[:,15],PatchYrCwtc0[:,15])
#plt.legend((expID[0:2]),loc=0,ncol=1) # title='PFT',
#plt.ylim((0.0,60))
plt.title('d', x=0.1, y=0.85,fontdict=font)
#plt.tick_params(labelbottom='off')
plt.xlabel('Year', fontdict=font)
plt.ylabel('Soil C (Kg C m$^{-2}$)', fontdict=font)


#%% Plotting yearly First cohort
xyear = RawData[0:yrmax,0].astype(np.float)
plt.figure(3) #
#plt.rc('axes', prop_cycle=(cycler('color', ['b','r', 'g',  'y', 'm'])+
#    cycler('linestyle', ['-', '--', ':', '-.',(0, (3, 1, 1, 1, 1, 1))])))
plt.clf()

plt.subplot(221)
plt.plot(xyear, CCYrVwtc0[:,8],CCYrCwtc0[:,8])
plt.legend((['Increasing WTC0', 'Constant WTC0']),loc=0,ncol=1) # title='PFT',
#plt.ylim((0.0,60))
plt.title('a', x=0.1, y=0.85,fontdict=font)
plt.tick_params(labelbottom='off')
#plt.xlabel('Year', fontdict=font)
plt.ylabel('Height(m)', fontdict=font)

plt.subplot(222)
plt.plot(xyear, CCYrVwtc0[:,6],CCYrCwtc0[:,6])
#plt.legend((expID[0:2]),loc=0,ncol=1) # title='PFT',
#plt.ylim((0.0,60))
plt.title('b', x=0.1, y=0.85,fontdict=font)
plt.tick_params(labelbottom='off')
#plt.xlabel('Year', fontdict=font)
plt.ylabel('dDBH (mm yr$^{-1}$)', fontdict=font)

plt.subplot(223)
plt.scatter(CCYrVwtc0[:,8], CCYrVwtc0[:,24], s=3)
#plt.legend((expID[0:2]),loc=0,ncol=1) # title='PFT',
#plt.ylim((0.0,60))
plt.title('c', x=0.1, y=0.85,fontdict=font)
#plt.tick_params(labelbottom='off')
plt.xlabel('Tree height (m)', fontdict=font)
plt.ylabel('Ktrunk\n(Kg H$_{2}$O/MPa/s)', fontdict=font)

plt.subplot(224)
plt.scatter(CCYrCwtc0[:,8], CCYrCwtc0[:,24], s=3,c='y')
#plt.legend((expID[0:2]),loc=0,ncol=1) # title='PFT',
#plt.ylim((0.0,60))
plt.title('d', x=0.1, y=0.85,fontdict=font)
#plt.tick_params(labelbottom='off')
plt.xlabel('Tree height (m)', fontdict=font)
plt.ylabel('Ktrunk\n(Kg H$_{2}$O/MPa/s)', fontdict=font)

#%% Plotting Hourly First cohort
# First cohort hourly dimensions
rows = len(CCHrV)
col  = len(CCHrV[1]) - 1
CCHrVwtc0 = CCHrV[0:rows,0:19].astype(np.float)
CCHrCwtc0 = CCHrC[0:rows,0:19].astype(np.float)
a_hour = np.arange(rows)
n1 = 20070
n2 = 20190

plt.figure(4) #
#plt.rc('axes', prop_cycle=(cycler('color', ['b','r', 'g',  'y', 'm'])+
#    cycler('linestyle', ['-', '--', ':', '-.',(0, (3, 1, 1, 1, 1, 1))])))
plt.clf()

plt.subplot(221)
plt.plot(a_hour[n1:n2], 1000*CCHrVwtc0[n1:n2,12]/CCHrVwtc0[n1:n2,9],a_hour[n1:n2], 1000*CCHrCwtc0[n1:n2,12]/CCHrCwtc0[n1:n2,9])
#plt.legend((['Increasing WTC0', 'Constant WTC0']),loc=0,ncol=1) # title='PFT',
plt.ylim((0.0,0.7))
plt.title('a', x=0.1, y=0.85,fontdict=font)
plt.tick_params(labelbottom='off')
#plt.xlabel('Year', fontdict=font)
plt.ylabel('GPP', fontdict=font)

plt.subplot(222)
plt.plot(a_hour[n1:n2], 1000*CCHrVwtc0[n1:n2,14]/CCHrVwtc0[n1:n2,9],a_hour[n1:n2], 1000*CCHrCwtc0[n1:n2,14]/CCHrCwtc0[n1:n2,9])
#plt.legend((expID[0:2]),loc=0,ncol=1) # title='PFT',
#plt.ylim((0.0,60))
plt.title('b', x=0.1, y=0.85,fontdict=font)
plt.tick_params(labelbottom='off')
#plt.xlabel('Year', fontdict=font)
plt.ylabel('Transp (mm yr$^{-1}$)', fontdict=font)

plt.subplot(223)
plt.plot(a_hour[n1:n2], CCHrVwtc0[n1:n2,15],a_hour[n1:n2], CCHrCwtc0[n1:n2,15])
#plt.legend((['Increasing WTC0', 'Constant WTC0']),loc=0,ncol=1) # title='PFT',
#plt.ylim((0.0,60))
plt.title('c', x=0.1, y=0.85,fontdict=font)
plt.tick_params(labelbottom='off')
plt.xlabel('Hour', fontdict=font)
plt.ylabel('$\psi$ $_{Leaf}$(MPa)', fontdict=font)

plt.subplot(224)
plt.plot(a_hour[n1:n2], CCHrVwtc0[n1:n2,16],a_hour[n1:n2], CCHrCwtc0[n1:n2,16])
#plt.legend((expID[0:2]),loc=0,ncol=1) # title='PFT',
#plt.ylim((0.0,60))
plt.title('d', x=0.1, y=0.85,fontdict=font)
plt.tick_params(labelbottom='off')
plt.xlabel('Hour', fontdict=font)
plt.ylabel('$\psi$ $_{SW}$(MPa)', fontdict=font)

#%% Plot Single PFT growth

fname1 = fpath+expID[3]+'/'+fname[5]+runID+'.csv'
f=open(fname1)
next(f) # skip headings
CCYrSingle=(np.array(list(csv.reader(f,delimiter=','))))
f.close()

# Cohort yearly dimensions, 5 PFTs
rows = len(CCYrSingle)
col  = len(CCYrSingle[1])
CCYrS = CCYrSingle[0:rows,0:col].astype(np.float)
n = int(rows/5)
CCYr5PFTs = CCYrS.reshape((n,5,col))
xyear  = CCYr5PFTs[:,4,0]
dDBH   = CCYr5PFTs[:,:,3]
Acrown = CCYr5PFTs[:,:,6]
GPP    = CCYr5PFTs[:,:,15]
Transp = CCYr5PFTs[:,:,17]
Atrunk = CCYr5PFTs[:,:,21]
Asap   = CCYr5PFTs[:,:,22]
Ktree  = CCYr5PFTs[:,:,23]
treeHU = CCYr5PFTs[:,:,24]
treeW0 = CCYr5PFTs[:,:,25]

plt.figure(5) #
#plt.rc('axes', prop_cycle=(cycler('color', ['b','r', 'g',  'y', 'm'])+
#    cycler('linestyle', ['-', '--', ':', '-.',(0, (3, 1, 1, 1, 1, 1))])))
plt.clf()

plt.subplot(321)
plt.plot(xyear, Transp/Acrown)
#plt.legend((expID[0:2]),loc=0,ncol=1) # title='PFT',
#plt.ylim((0.0,60))
plt.title('a', x=0.1, y=0.85,fontdict=font)
plt.tick_params(labelbottom='off')
#plt.xticks(np.arange(0, 400, 100))
#plt.xlabel('Year', fontdict=font)
plt.ylabel('Transp (mm yr$^{-1}$)', fontdict=font)

plt.subplot(322)
plt.plot(xyear, dDBH)
#plt.legend((expID[0:2]),loc=0,ncol=1) # title='PFT',
#plt.ylim((0.0,60))
plt.title('b', x=0.1, y=0.85,fontdict=font)
plt.tick_params(labelbottom='off')
#plt.xticks(np.arange(0, 400, 100))
#plt.xlabel('Year', fontdict=font)
plt.ylabel('dDBH (mm yr$^{-1}$)', fontdict=font)

plt.subplot(323)
#plt.plot(xyear, Ktree) #100*Ktree/Acrown
plt.plot(xyear, 100*Ktree/Acrown)
#plt.legend((expID[0:2]),loc=0,ncol=1) # title='PFT',
#plt.ylim((0.0,60))
plt.title('c', x=0.1, y=0.85,fontdict=font)
plt.tick_params(labelbottom='off')
#plt.xticks(np.arange(0, 400, 100))
#plt.xlabel('Year', fontdict=font)
plt.ylabel('Ktrunk (m$^{3}$/MPa/s)', fontdict=font)

plt.subplot(324)
plt.plot(xyear, Asap/Acrown*10000/3.5)
#plt.ylim((0.0,60))
plt.title('d', x=0.1, y=0.85,fontdict=font)
plt.tick_params(labelbottom='off')
#plt.ylabel('GPP\n(Kg C m$^{-2}$ yr$^{-1}$)', fontdict=font)
#plt.xticks(np.arange(0, 400, 100))
#plt.xlabel('Year', fontdict=font)
plt.ylabel('Hv (cm$^{2}$/m$^{2}$)', fontdict=font)

plt.subplot(325)
plt.plot(xyear, Asap/Atrunk)
plt.ylim((0.0,0.8))
plt.title('e', x=0.1, y=0.85,fontdict=font)
plt.tick_params(labelbottom='off')
#plt.ylabel('GPP\n(Kg C m$^{-2}$ yr$^{-1}$)', fontdict=font)
#plt.xticks(np.arange(0, 400, 100))
plt.xlabel('Year', fontdict=font)
plt.ylabel('$r_{SW/Trunk}$', fontdict=font)

plt.subplot(326)
plt.plot(xyear, treeHU/treeW0)
plt.legend((PFTID),loc=0,ncol=2) # title='PFT',
plt.ylim((0.0,0.95))
plt.title('f', x=0.1, y=0.85,fontdict=font)
#plt.tick_params(labelbottom='off')
#plt.xticks(np.arange(0, 400, 100))
plt.xlabel('Year', fontdict=font)
plt.ylabel('$r_{HU/W0}$', fontdict=font)
