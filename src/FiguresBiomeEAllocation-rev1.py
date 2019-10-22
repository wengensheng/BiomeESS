#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 21:06:50 2017

@author: eweng
"""
import numpy as np
import csv
import matplotlib.pyplot as plt

font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 12,
        }

# "Tests" file, for monoculture runs, filename1
#0	  1	2	3	4	5	6	
#year,BA1 BA2         	BA3         	BA4         	BA5         	BA6         	
#7	8	9	10	11	12	
#GPP1  GPP2  GPP3  GPP4  GPP5  GPP6
#13	14	15	16	17	18	
#NPP1  NPP2  NPP3  NPP4  NPP5  NPP6        	
#19	20	21	22	23	24	
#NPPlf1      	NPPlf2      	NPPlf3      	NPPlf4      	NPPlf5      	NPPlf6
#25	26	27	28	29	30	
#NPPfr1      	NPPfr2      	NPPfr3      	NPPfr4      	NPPfr5      	NPPfr6
#31	32	33	34	35	36	
#NPPwd1      	NPPwd2      	NPPwd3      	NPPwd4      	NPPwd5      	NPPwd6
#37	38	39	40	41	42	
#NPPsd1      	NPPsd2      	NPPsd3      	NPPsd4      	NPPsd5      	NPPsd6
#43	44	45	46	47	48	
#NSC1        	NSC2        	NSC3        	NSC4        	NSC5        	NSC6
#49	50	51	52	53	54	
#NSN1        	NSN2        	NSN3        	NSN4        	NSN5        	NSN6
#55	56	57	58	59	60	
#CAI-1       	CAI-2       	CAI-3       	CAI-4       	CAI-5       	CAI-6
#61	62	63	64	65	66	
#HTcri1      	HTcri2      	HTcri3      	HTcri4      	HTcri5      	HTcri6
#67	68	69	70	71	72	
#PlantC1     	PlantC2     	PlantC3     	PlantC4     	PlantC5     	PlantC6
#73	74	75	76	77	78	
#SoilC1      	SoilC2      	SoilC3      	SoilC4      	SoilC5      	SoilC6
#79	80	81	82	83	84	
#PlantN1     	PlantN2     	PlantN3     	PlantN4     	PlantN5     	PlantN6
#85	86	87	88	89	90	
#SoilN1      	SoilN2      	SoilN3      	SoilN4      	SoilN5      	SoilN6
#91	92	93	94	95	96
#Nmnr1       	Nmnr2       	Nmnr3       	Nmnr4       	Nmnr5       	Nmnr6       

# 'Sum' files, for multiculture runs, filename2
#0	  1	2	3	4	5	6	
#year,BA1,BA2,BA3,BA4,BA5,BA6         	
#7	8	9	10	11	12	
#NPPlf1,NPPlf2,NPPlf3,NPPlf4,NPPlf5 NPPlf6      	
#13	14	15	16	17	18 	
#NPPfr1      	NPPfr2      	NPPfr3      	NPPfr4      	NPPfr5      	NPPfr6
#19	20	21	22	23	24
#NPPwd1      	NPPwd2      	NPPwd3      	NPPwd4      	NPPwd5      	NPPwd6
#25	26	27	28	29	30	
#NPPsd1      	NPPsd2      	NPPsd3      	NPPsd4      	NPPsd5      	NPPsd6 
#31	32	33	34	35	36	  	
#dDBH1       	dDBH2       	dDBH3       	dDBH4       	dDBH5       	dDBH6  
#37	38	39	40 41	42	
#N_PFT1      	N_PFT2      	N_PFT3      	N_PFT4      	N_PFT5      	N_PFT6 
#43	44	45	46	47	48	     	
#cLAI1       	cLAI2       	cLAI3       	cLAI4       	cLAI5       	cLAI6
#49	50	51	52	53	54	55	56	57	58	
#DBH0.1, DBH0.2, DBH0.3, DBH0.4, DBH0.5, DBH0.6, DBH0.7, DBH0.8, DBH0.9, DBH1.0
#59	60	61	62	63		      	     
#siteNPP	   siteNPPLF	   siteNPPFR	   siteNPPWD	   siteNPPsd	
# 64	   65  66	67	       
#totBA	       critH	         CAI	         LAI	 
#	68	69     70		        
#   GPP	Rauto  Rh	
# 71	    72   73	74	75		       
#rain	   SiolWater	      Transp	        Evap	      Runoff	
#76	    77   78	79	80		      
#plantC	    soilC	      plantN	       soilN	        totN	
#81	82	83   84	85	86         
#NSC	       SeedC	       leafC	       rootC	    SapwoodC	       WoodC	
#	87	88    89	90	91	92        
#NSN	       SeedN	       leafN	       rootN	    SapwoodN	       WoodN	
#93	94  95	96	97	98	99	       
#McrbC, fastSOM, SlowSOM, McrbN, fastSoilN, slowSoilN	, mineralN, 
# 100	101     102
#N_fxed  N_uptk  N_yrMin

datapath1 = 'model/output/rerun0521/SC0.1/sum0523/' #'model/output/forPaper/'
#datapath2 = 'model/output/rerun0820/test10/' # 'csv/forPaper/'
datapath2 = 'model/output/rerun0521/SC0.1/sum0523/' #'model/output/rerun1120/'

#%% Ambient CO2
runtag = '_aCO2_' #  '_eCO2_' # '_fN2FixA_'

# Single PFT runs

filename1 = ['Tests'+runtag+'SC04_One.csv', \
             'Tests'+runtag+'SC06_One.csv', \
             'Tests'+runtag+'SC08_One.csv', \
             'Tests'+runtag+'SC10_One.csv', \
             'Tests'+runtag+'SC12_One.csv', \
             'Tests'+runtag+'SC14_One.csv', \
             'Tests'+runtag+'SC16_One.csv', \
             'Tests'+runtag+'SC18_One.csv', \
             ]

# All PFTs runs

filename2 = ['Sum_SC04'+runtag+'All.csv', \
             'Sum_SC06'+runtag+'All.csv', \
             'Sum_SC08'+runtag+'All.csv', \
             'Sum_SC10'+runtag+'All.csv', \
             'Sum_SC12'+runtag+'All.csv', \
             'Sum_SC14'+runtag+'All.csv', \
             'Sum_SC16'+runtag+'All.csv', \
             'Sum_SC18'+runtag+'All.csv', \
             ]
# Data used for making figures
Site_avg =np.zeros((8,102)) # for mono runs

# for poly runs, 8 pfts
BAall = np.zeros((8,1800,8)) # 8 N levels, 1800 years, 8 PFTs
LAIsite = np.zeros((1800,8))
BMsite = np.zeros((1800,8))
SNsite = np.zeros((1800,8))
HTsite = np.zeros((1800,8))
GPPsite = np.zeros((1800,8))
NPPsite = np.zeros((1800,8))
NPLsite = np.zeros((1800,8))
NPRsite = np.zeros((1800,8))
NPWsite = np.zeros((1800,8))
minNsite = np.zeros((1800,8))

# 7: the first 6 for the mono R/L runs, the 7th for poly run
HT_avg  = np.zeros((8,7))
BM_avg  = np.zeros((8,7))
HT_avg  = np.zeros((8,7))
GPP_avg = np.zeros((8,7))
NPP_avg = np.zeros((8,7))
NPL_avg = np.zeros((8,7))
NPR_avg = np.zeros((8,7))
NPW_avg = np.zeros((8,7))
Nm_avg = np.zeros((8,7))
LAI_avg = np.zeros((8,7))

#Nlevels = np.array([130,188,245,302,359,416])
Nlevels = np.array([114.5,177,239.5,302,364.5,427,489.5,552])
RLlevels = np.array([1,2,3,4,5,6])
#datapath2 = 'model/output/rerun1120/'
for i in range(8): # nitrogen levels
    f=open(datapath2 +filename1[i]) # mono runs
    next(f) # skip headings
    sitedata=list(csv.reader(f,delimiter=','))
    f.close()
    
    ecodatayr=np.array(sitedata)
    xyear = ecodatayr[0:1800,0].astype(np.float)
    Site_all = ecodatayr[0:1800,1:103].astype(np.float)
    Site_avg[i,:] = np.mean(Site_all[1400:1800,:],axis=0)

#datapath2 = 'model/output/rerun1120/' #FixedRL8PFTs/'
for i in range(8):
    f=open(datapath2 +filename2[i])
    next(f) # skip headings
    yearlydata=list(csv.reader(f,delimiter=','))
    f.close()
    ecodatayr=np.array(yearlydata)
    Sum_all = ecodatayr[0:1800,1:102].astype(np.float)
    xyear = ecodatayr[0:1800,0].astype(np.float)
    BAall[i,:,:] = ecodatayr[0:1800,1:9].astype(np.float)
    LAIsite[:,i] = ecodatayr[0:1800,83].astype(np.float)
    GPPsite[:,i] = ecodatayr[0:1800,84].astype(np.float) # site NPP
    NPPsite[:,i] = ecodatayr[0:1800,75].astype(np.float) # site NPP
    NPLsite[:,i] = ecodatayr[0:1800,76].astype(np.float)
    NPRsite[:,i] = ecodatayr[0:1800,78].astype(np.float) # site root NPP
    NPWsite[:,i] = ecodatayr[0:1800,77].astype(np.float) # site woodNPP
    BMsite[:,i] = ecodatayr[0:1800,92].astype(np.float) # Plant biomass
    HTsite[:,i] = ecodatayr[0:1800,81].astype(np.float) # Critical height
    minNsite[:,i] = ecodatayr[0:1800,115].astype(np.float) # mineral N
# multi-year mean for multiculture runs
siteBMavg = np.mean(BMsite[1400:1800,:],axis=0)
siteHTavg = np.mean(HTsite[1200:1800,:],axis=0)
siteLAIavg = np.mean(LAIsite[1400:1800,:],axis=0)
siteGPPavg = np.mean(GPPsite[1400:1800,:],axis=0)
siteNPPavg = np.mean(NPPsite[1400:1800,:],axis=0)
siteNPLavg = np.mean(NPLsite[1400:1800,:],axis=0)
siteNPRavg = np.mean(NPRsite[1400:1800,:],axis=0)
siteNPWavg = np.mean(NPWsite[1400:1800,:],axis=0)
siteNmavg = np.mean(minNsite[1400:1800,:],axis=0)

# Monoculture results, 6 PFTs in total
HT_avg[:,0:6]  = Site_avg[:,60:66]
BM_avg[:,0:6]  = Site_avg[:,66:72]
GPP_avg[:,0:6] = Site_avg[:,6:12]/1000.
NPP_avg[:,0:6] = Site_avg[:,12:18]/1000.
NPL_avg[:,0:6] = Site_avg[:,18:24]/1000.
NPR_avg[:,0:6] = Site_avg[:,24:30]/1000.
NPW_avg[:,0:6] = Site_avg[:,30:36]/1000.
Nm_avg[:,0:6]  = Site_avg[:,90:96]
LAI_avg[:,0:6] = Site_avg[:,96:102]

#Multiculture runs, 7th column
BM_avg[:,6]  = np.transpose(siteBMavg)
HT_avg[:,6]  = np.transpose(siteHTavg)
GPP_avg[:,6] = np.transpose(siteGPPavg)
NPP_avg[:,6] = np.transpose(siteNPPavg)/1000.
NPL_avg[:,6] = np.transpose(siteNPLavg)/1000.
NPR_avg[:,6] = np.transpose(siteNPRavg)/1000.
NPW_avg[:,6] = np.transpose(siteNPWavg)/1000.
Nm_avg[:,6] = np.transpose(siteNmavg)
LAI_avg[:,6] = np.transpose(siteLAIavg)
 
np.savetxt('MonocultureAverage_aCO2.csv', np.transpose(Site_avg), delimiter=',') 
np.savetxt('BiomassAverage_aCO2.csv', BM_avg, delimiter=',')
GPPaCO2 = np.copy(GPP_avg)
NPPaCO2 = np.copy(NPP_avg)
BAaCO2 = np.copy(BAall)
BMaCO2 = np.copy(BM_avg)
HTaCO2 = np.copy(HT_avg)
NPLaCO2 = np.copy(NPL_avg)
NPRaCO2 = np.copy(NPR_avg)
NPWaCO2 = np.copy(NPW_avg)
NmaCO2 = np.copy(Nm_avg)
LAIaCO2 = np.copy(LAI_avg)
#%% Plot
xmax = 1800
ymax = 50



#%% elevated CO2
runtag = '_eCO2_' 
# Single PFT runs

filename1 = ['Tests'+runtag+'SC04_One.csv', \
             'Tests'+runtag+'SC06_One.csv', \
             'Tests'+runtag+'SC08_One.csv', \
             'Tests'+runtag+'SC10_One.csv', \
             'Tests'+runtag+'SC12_One.csv', \
             'Tests'+runtag+'SC14_One.csv', \
             'Tests'+runtag+'SC16_One.csv', \
             'Tests'+runtag+'SC18_One.csv', \
             ]

# All PFTs runs

filename2 = ['Sum_SC04'+runtag+'All.csv', \
             'Sum_SC06'+runtag+'All.csv', \
             'Sum_SC08'+runtag+'All.csv', \
             'Sum_SC10'+runtag+'All.csv', \
             'Sum_SC12'+runtag+'All.csv', \
             'Sum_SC14'+runtag+'All.csv', \
             'Sum_SC16'+runtag+'All.csv', \
             'Sum_SC18'+runtag+'All.csv', \
             ]
BAall = np.zeros((8,1800,8))
#datapath2 = 'model/output/rerun1120/'
for i in range(8):
    f=open(datapath2 +filename1[i])
    next(f) # skip headings
    sitedata=list(csv.reader(f,delimiter=','))
    f.close()
    ecodatayr=np.array(sitedata)
    xyear = ecodatayr[0:1800,0].astype(np.float)
    Site_all = ecodatayr[0:1800,1:103].astype(np.float)
    Site_avg[i,:] = np.mean(Site_all[1400:1800,:],axis=0)

#datapath2 = 'model/output/rerun1120/' #FixedRL8PFTs/'
for i in range(8):
    f=open(datapath2+filename2[i])
    next(f) # skip headings
    yearlydata=list(csv.reader(f,delimiter=','))
    f.close()
    ecodatayr=np.array(yearlydata)
    xyear = ecodatayr[0:1800,0].astype(np.float)
    BAall[i,:,:] = ecodatayr[0:1800,1:9].astype(np.float)
    LAIsite[:,i] = ecodatayr[0:1800,83].astype(np.float)
    GPPsite[:,i] = ecodatayr[0:1800,84].astype(np.float) # site NPP
    NPPsite[:,i] = ecodatayr[0:1800,75].astype(np.float) # site NPP
    NPLsite[:,i] = ecodatayr[0:1800,76].astype(np.float)
    NPRsite[:,i] = ecodatayr[0:1800,78].astype(np.float) # site root NPP
    NPWsite[:,i] = ecodatayr[0:1800,77].astype(np.float) # site woodNPP
    BMsite[:,i] = ecodatayr[0:1800,92].astype(np.float) # Plant biomass
    HTsite[:,i] = ecodatayr[0:1800,81].astype(np.float) # Critical height
    minNsite[:,i] = ecodatayr[0:1800,115].astype(np.float) # mineral N
# multi-year mean for multiculture runs
siteBMavg = np.mean(BMsite[1400:1800,:],axis=0)
siteHTavg = np.mean(HTsite[1200:1800,:],axis=0)
siteLAIavg = np.mean(LAIsite[1400:1800,:],axis=0)
siteGPPavg = np.mean(GPPsite[1400:1800,:],axis=0)
siteNPPavg = np.mean(NPPsite[1400:1800,:],axis=0)
siteNPLavg = np.mean(NPLsite[1400:1800,:],axis=0)
siteNPRavg = np.mean(NPRsite[1400:1800,:],axis=0)
siteNPWavg = np.mean(NPWsite[1400:1800,:],axis=0)    
siteNmavg = np.mean(minNsite[1400:1800,:],axis=0)

#%%
# 7: the first 6 for the mono runs, the 7th for poly run
BM_avg  = np.zeros((8,7))
HT_avg  = np.zeros((8,7))
GPP_avg = np.zeros((8,7))
NPP_avg = np.zeros((8,7))
NPL_avg = np.zeros((8,7))
NPR_avg = np.zeros((8,7))
NPW_avg = np.zeros((8,7))
Nm_avg = np.zeros((8,7))
LAI_avg = np.zeros((8,7))
# Monoculture results, 6 PFTs in total
HT_avg[:,0:6] = Site_avg[:,60:66]
BM_avg[:,0:6] = Site_avg[:,66:72]
GPP_avg[:,0:6] = Site_avg[:,6:12]/1000.
NPP_avg[:,0:6] = Site_avg[:,12:18]/1000.
NPL_avg[:,0:6] = Site_avg[:,18:24]/1000.
NPR_avg[:,0:6] = Site_avg[:,24:30]/1000.
NPW_avg[:,0:6] = Site_avg[:,30:36]/1000.
Nm_avg[:,0:6]  = Site_avg[:,90:96]
LAI_avg[:,0:6] = Site_avg[:,96:102]
#Multiculture runs, 7th column
BM_avg[:,6]  = np.transpose(siteBMavg)
HT_avg[:,6]  = np.transpose(siteHTavg)
GPP_avg[:,6] = np.transpose(siteGPPavg)
NPP_avg[:,6] = np.transpose(siteNPPavg)/1000.
NPL_avg[:,6] = np.transpose(siteNPLavg)/1000.
NPR_avg[:,6] = np.transpose(siteNPRavg)/1000.
NPW_avg[:,6] = np.transpose(siteNPWavg)/1000.
Nm_avg[:,6] = np.transpose(siteNmavg)
LAI_avg[:,6] = np.transpose(siteLAIavg)

np.savetxt('MonocultureAverage_eCO2.csv', np.transpose(Site_avg), delimiter=',') 
np.savetxt('BiomassAverage_eCO2.csv', BM_avg, delimiter=',') 
BAeCO2  = np.copy(BAall)
BMeCO2  = np.copy(BM_avg)
HTeCO2  = np.copy(HT_avg)
GPPeCO2 = np.copy(GPP_avg)
NPPeCO2 = np.copy(NPP_avg)
NPLeCO2 = np.copy(NPL_avg)
NPReCO2 = np.copy(NPR_avg)
NPWeCO2 = np.copy(NPW_avg)
NmeCO2  = np.copy(Nm_avg)
LAIeCO2 = np.copy(LAI_avg)

#%% Plot
# Mono GPP, NPP, Allocation, Plant biomass #Ray's suggestion
plt.figure(1)
plt.subplot(321)
plt.plot(Nlevels,(GPPaCO2[:,0:6]),'o--',mfc='none')
plt.title('(a)', x=0.1, y=0.85,fontdict=font)
#plt.xlim((0,6))
plt.ylim((0,1.7))
plt.ylabel('GPP\n(kgC $m^{-2} yr^{-1}$)', fontdict=font)

plt.subplot(322)
plt.plot(Nlevels,(NPPaCO2[:,0:6]),'o--',mfc='none')
plt.title('(b)', x=0.1, y=0.85,fontdict=font)
plt.ylim((0,1.7))
plt.ylabel('NPP\n(kgC $m^{-2} yr^{-1}$)', fontdict=font)
#plt.legend( ('130','188','245',\
#             '302','359','416'),\
#             title="Total N ($gN/m^{2}$)",loc=0, ncol=2)
plt.legend( ('$\phi_{RL}$: 1', '$\phi_{RL}$: 2', '$\phi_{RL}$: 3', \
             '$\phi_{RL}$: 4', '$\phi_{RL}$: 5', '$\phi_{RL}$: 6'),\
              loc=0,ncol=2)
plt.subplot(323)
plt.plot(Nlevels,(NPRaCO2[:,0:6]/NPPaCO2[:,0:6]),'o--',mfc='none')
plt.title('(c)', x=0.1, y=0.85,fontdict=font)
plt.ylim((0,1.0))
plt.ylabel('$f_{NPP,Root}$', fontdict=font)

plt.subplot(324)
plt.plot(Nlevels,(NPWaCO2[:,0:6]/NPPaCO2[:,0:6]),'o--',mfc='none')
plt.title('(d)', x=0.1, y=0.85,fontdict=font)
plt.ylim((0,1.0))
plt.ylabel('$f_{NPP,Wood}$', fontdict=font)
#plt.xlabel('Root/leaf area ratio ($\phi_{RL}$)', fontdict=font)

plt.subplot(325)
plt.plot(Nlevels,(NPLaCO2[:,0:6]/NPPaCO2[:,0:6]),'o--',mfc='none')
plt.title('(e)', x=0.1, y=0.85,fontdict=font)
plt.ylim((0.1,0.18))
plt.ylabel('$f_{NPP,Leaf}$', fontdict=font)
plt.xlabel('Total N ($gN/m^{2}$)', fontdict=font)

plt.subplot(326)
plt.plot(Nlevels,(BMaCO2[:,0:6]),'o--',mfc='none')
plt.title('(f)', x=0.1, y=0.85,fontdict=font)
plt.ylim((0,30))
plt.ylabel('Plant biomass\n(kgC $m^{-2}$)', fontdict=font)
plt.xlabel('Total N ($gN/m^{2}$)', fontdict=font)

# Monoculture results of eCO2
plt.figure(2)
plt.subplot(321)
plt.plot(Nlevels,(GPPeCO2[:,0:6]),'o--',mfc='none')
plt.title('(a)', x=0.1, y=0.85,fontdict=font)
#plt.xlim((0,6))
plt.ylim((0,2.5))
plt.ylabel('GPP\n(kgC $m^{-2} yr^{-1}$)', fontdict=font)

plt.subplot(322)
plt.plot(Nlevels,(NPPeCO2[:,0:6]),'o--',mfc='none')
plt.title('(b)', x=0.1, y=0.85,fontdict=font)
plt.ylim((0,2.5))
plt.ylabel('NPP\n(kgC $m^{-2} yr^{-1}$)', fontdict=font)
#plt.legend( ('130','188','245',\
#             '302','359','416'),\
#             title="Total N ($gN/m^{2}$)",loc=0, ncol=2)
plt.legend( ('$\phi_{RL}$: 1', '$\phi_{RL}$: 2', '$\phi_{RL}$: 3', \
             '$\phi_{RL}$: 4', '$\phi_{RL}$: 5', '$\phi_{RL}$: 6'),\
              loc=0,ncol=2)
plt.subplot(323)
plt.plot(Nlevels,(NPReCO2[:,0:6]/NPPeCO2[:,0:6]),'o--',mfc='none')
plt.title('(c)', x=0.1, y=0.85,fontdict=font)
plt.ylim((0,1.0))
plt.ylabel('$f_{NPP,Root}$', fontdict=font)

plt.subplot(324)
plt.plot(Nlevels,(NPWeCO2[:,0:6]/NPPeCO2[:,0:6]),'o--',mfc='none')
plt.title('(d)', x=0.1, y=0.85,fontdict=font)
plt.ylim((0,1.0))
plt.ylabel('$f_{NPP,Wood}$', fontdict=font)
#plt.xlabel('Root/leaf area ratio ($\phi_{RL}$)', fontdict=font)

plt.subplot(325)
plt.plot(Nlevels,(NPLeCO2[:,0:6]/NPPeCO2[:,0:6]),'o--',mfc='none')
plt.title('(e)', x=0.1, y=0.85,fontdict=font)
plt.ylim((0.08,0.16))
plt.ylabel('$f_{NPP,Leaf}$', fontdict=font)
plt.xlabel('Total N ($gN/m^{2}$)', fontdict=font)

plt.subplot(326)
plt.plot(Nlevels,(BMeCO2[:,0:6]),'o--',mfc='none')
plt.title('(f)', x=0.1, y=0.85,fontdict=font)
plt.ylim((0,50))
plt.ylabel('Plant biomass\n(kgC $m^{-2}$)', fontdict=font)
plt.xlabel('Total N ($gN/m^{2}$)', fontdict=font)
#%% Scuccession
#aCO2 succession
#[114.5,177,239.5,302,364.5,427,489.5,552]
plt.figure(3)
plt.subplot(321)
plt.plot(xyear, BAaCO2[0,:,0:6])
plt.xlim((0,xmax))
plt.ylim((0,30))
plt.title('114.5 $gN/m^{2}$', x=0.5, y=0.85, fontdict=font)
#plt.xlabel('Year', fontdict=font)
plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)
#plt.show()

plt.subplot(322)
plt.plot(xyear, BAaCO2[1,:,0:6])
plt.xlim((0,xmax))
plt.ylim((0,30))
plt.legend( ('$\phi_{RL}$: 1', '$\phi_{RL}$: 2', '$\phi_{RL}$: 3', \
             '$\phi_{RL}$: 4', '$\phi_{RL}$: 5', '$\phi_{RL}$: 6', \
             '$\phi_{RL}$: 7', '$\phi_{RL}$: 8'), \
             ncol=1, loc='center left', bbox_to_anchor=(1, 0.5) )
plt.title('177 $gN/m^{2}$',x=0.5, y=0.85,  fontdict=font)
#plt.xlabel('Year', fontdict=font)
#plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)

plt.subplot(323)
plt.plot(xyear, BAaCO2[2,:,0:6])
plt.xlim((0,xmax))
plt.ylim((0,30))
#plt.legend( ('R/L: 1', 'R/L: 2', 'R/L: 3','R/L: 4', 'R/L: 5') )
plt.title('239.5 $gN/m^{2}$', x=0.5, y=0.85,fontdict=font)
#plt.xlabel('Year', fontdict=font)
plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)

plt.subplot(324)
plt.plot(xyear, BAaCO2[3,:,0:6])
plt.xlim((0,xmax))
plt.ylim((0,ymax))
#plt.legend( ('R/L: 1', 'R/L: 2', 'R/L: 3','R/L: 4', 'R/L: 5') )
plt.title('302 $gN/m^{2}$', x=0.5, y=0.85,fontdict=font)
#plt.xlabel('Year', fontdict=font)
#plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)

plt.subplot(325)
plt.plot(xyear, BAaCO2[4,:,0:6])
plt.xlim((0,xmax))
plt.ylim((0,ymax))
#plt.legend( ('R/L: 1', 'R/L: 2', 'R/L: 3','R/L: 4', 'R/L: 5') )
plt.title('364.5 $gN/m^{2}$', x=0.5, y=0.85,fontdict=font)
plt.xlabel('Year', fontdict=font)
plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)

plt.subplot(326)
plt.plot(xyear, BAaCO2[5,:,0:6])
plt.xlim((0,xmax))
plt.ylim((0,ymax))
#plt.legend( ('R/L: 1', 'R/L: 2', 'R/L: 3','R/L: 4', 'R/L: 5') )
plt.title('427 $gN/m^{2}$', x=0.5, y=0.85,fontdict=font)
plt.xlabel('Year', fontdict=font)
#plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)

# Succession eCO2
plt.figure(4)
plt.subplot(321)
plt.plot(xyear, BAeCO2[0,:,0:8])
plt.xlim((0,xmax))
plt.ylim((0,ymax))
plt.title('114.5 $gN/m^{2}$', x=0.5, y=0.85, fontdict=font)
#plt.xlabel('Year', fontdict=font)
plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)
#plt.show()

plt.subplot(322)
plt.plot(xyear, BAeCO2[1,:,0:8])
plt.xlim((0,xmax))
plt.ylim((0,ymax))
plt.legend( ('$\phi_{RL}$: 1', '$\phi_{RL}$: 2', '$\phi_{RL}$: 3', \
             '$\phi_{RL}$: 4', '$\phi_{RL}$: 5', '$\phi_{RL}$: 6', \
             '$\phi_{RL}$: 7', '$\phi_{RL}$: 8'), \
             ncol=1, loc='center left', bbox_to_anchor=(1, 0.5) )
plt.title('177 $gN/m^{2}$',x=0.5, y=0.85,  fontdict=font)
#plt.xlabel('Year', fontdict=font)
#plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)

plt.subplot(323)
plt.plot(xyear, BAeCO2[2,:,0:8])
plt.xlim((0,xmax))
plt.ylim((0,ymax))
#plt.legend( ('R/L: 1', 'R/L: 2', 'R/L: 3','R/L: 4', 'R/L: 5') )
plt.title('239.5 $gN/m^{2}$', x=0.5, y=0.85,fontdict=font)
#plt.xlabel('Year', fontdict=font)
plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)

plt.subplot(324)
plt.plot(xyear, BAeCO2[3,:,0:8])
plt.xlim((0,xmax))
plt.ylim((0,ymax))
plt.title('302 $gN/m^{2}$', x=0.5, y=0.85,fontdict=font)
#plt.xlabel('Year', fontdict=font)
#plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)

plt.subplot(325)
plt.plot(xyear, BAeCO2[4,:,0:8])
plt.xlim((0,xmax))
plt.ylim((0,ymax))
#plt.legend( ('R/L: 1', 'R/L: 2', 'R/L: 3','R/L: 4', 'R/L: 5') )
plt.title('364.5 $gN/m^{2}$', x=0.5, y=0.85,fontdict=font)
plt.xlabel('Year', fontdict=font)
plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)

plt.subplot(326)
plt.plot(xyear, BAeCO2[5,:,0:8])
plt.xlim((0,xmax))
plt.ylim((0,ymax))
#plt.legend( ('R/L: 1', 'R/L: 2', 'R/L: 3','R/L: 4', 'R/L: 5') )
plt.title('427 $gN/m^{2}$', x=0.5, y=0.85,fontdict=font)
plt.xlabel('Year', fontdict=font)
#plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)

# Compare succession: aCO2 vs. eCO2
plt.figure(41)
plt.subplot(321)
plt.plot(xyear, BAaCO2[1,:,0:6])
plt.tick_params(labelbottom='off')
plt.xlim((0,xmax))
plt.ylim((0,ymax))
#plt.legend( ('$\phi_{RL}$: 1', '$\phi_{RL}$: 2', '$\phi_{RL}$: 3', \
#             '$\phi_{RL}$: 4', '$\phi_{RL}$: 5') )
plt.title('(a) 177 $gN/m^{2}$, aCO2', x=0.4, y=0.8, fontdict=font)
#plt.xlabel('Year', fontdict=font)
#plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)
#plt.show()

plt.subplot(323)
plt.plot(xyear, BAaCO2[3,:,0:6])
plt.tick_params(labelbottom='off')
plt.xlim((0,xmax))
plt.ylim((0,ymax))
#plt.legend( ('R/L: 1', 'R/L: 2', 'R/L: 3','R/L: 4', 'R/L: 5') )
plt.title('(c) 302 $gN/m^{2}$, aCO2', x=0.4, y=0.8,fontdict=font)
#plt.xlabel('Year', fontdict=font)
plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)

plt.subplot(325)
plt.plot(xyear, BAaCO2[5,:,0:6])
plt.xlim((0,xmax))
plt.ylim((0,ymax))
#plt.legend( ('R/L: 1', 'R/L: 2', 'R/L: 3','R/L: 4', 'R/L: 5') )
plt.title('(e) 427 $gN/m^{2}$, aCO2', x=0.4, y=0.8,fontdict=font)
plt.xlabel('Year', fontdict=font)
#plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)

plt.subplot(322)
plt.plot(xyear, BAeCO2[1,:,0:8])
plt.tick_params(labelbottom='off', labelleft='off')
plt.xlim((0,xmax))
plt.ylim((0,ymax))
plt.title('(b) 177 $gN/m^{2}$, eCO2', x=0.4, y=0.8, fontdict=font)
#plt.xlabel('Year', fontdict=font)
#plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)
#plt.show()

plt.subplot(324)
plt.plot(xyear, BAeCO2[3,:,0:8])
plt.tick_params(labelbottom='off', labelleft='off')
plt.xlim((0,xmax))
plt.ylim((0,ymax))
plt.legend( ('$\phi_{RL}$: 1', '$\phi_{RL}$: 2', '$\phi_{RL}$: 3', \
             '$\phi_{RL}$: 4', '$\phi_{RL}$: 5', '$\phi_{RL}$: 6', \
             '$\phi_{RL}$: 7', '$\phi_{RL}$: 8'), \
             ncol=1, loc='center left', bbox_to_anchor=(1, 0.5) )
plt.title('(d) 302 $gN/m^{2}$, eCO2', x=0.4, y=0.8,fontdict=font)
#plt.xlabel('Year', fontdict=font)
#plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)

plt.subplot(326)
plt.plot(xyear, BAeCO2[5,:,0:8])
plt.tick_params( labelleft='off')
plt.xlim((0,xmax))
plt.ylim((0,ymax))
plt.title('(f) 427 $gN/m^{2}$, eCO2', x=0.4, y=0.8,fontdict=font)
plt.xlabel('Year', fontdict=font)
#plt.ylabel('Basal area ($m^{2} ha^{-1}$)', fontdict=font)

#%% Mean results compare
plt.figure(51) # GPP and NPP of poly-runs

plt.subplot(311)
plt.plot(Nlevels, GPPaCO2[:,6],'bo-',label='Poly, aCO$_{2}$')
plt.plot(Nlevels, GPPeCO2[:,6],'rD-',label='Poly, eCO$_{2}$')
plt.plot(Nlevels, GPPaCO2[:,3],'bo--',mfc='none',label='$\phi_{RL}$=4, aCO$_{2}$')
plt.plot(Nlevels, GPPeCO2[:,3],'rD--',mfc='none',label='$\phi_{RL}$=4, eCO$_{2}$')
plt.legend(ncol=1, loc='lower right')
#plt.xlim((100,450))
plt.yticks(np.arange(0, 2.5, 0.5))
plt.title('(a)', x=0.1, y=0.85,fontdict=font)
#plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('GPP (kgC $m^{-2}$)', fontdict=font)

plt.subplot(312)
plt.plot(Nlevels, NPPaCO2[:,6],'bo-')
plt.plot(Nlevels, NPPeCO2[:,6],'rD-')
plt.plot(Nlevels, NPPaCO2[:,3],'bo--',mfc='none')
plt.plot(Nlevels, NPPeCO2[:,3],'rD--',mfc='none')
plt.legend()
#plt.xlim((100,450))
plt.yticks(np.arange(0, 1.5, 0.5))
plt.title('(b)', x=0.1, y=0.85,fontdict=font)
#plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
#plt.ylabel('$f_{NPP,leaf}$', fontdict=font)
plt.ylabel('NPP (kgC $m^{-2}$)', fontdict=font)

plt.subplot(313)
plt.plot(Nlevels, NPPaCO2[:,6]/GPPaCO2[:,6],'bo-')
plt.plot(Nlevels, NPPeCO2[:,6]/GPPeCO2[:,6],'rD-')
plt.plot(Nlevels, NPPaCO2[:,3]/GPPaCO2[:,3],'bo--',mfc='none')
plt.plot(Nlevels, NPPeCO2[:,3]/GPPeCO2[:,3],'rD--',mfc='none')
plt.legend()
#plt.xlim((100,450))
plt.ylim((0.5,0.6))
plt.title('(c)', x=0.1, y=0.85,fontdict=font)
plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('NPP/GPP', fontdict=font)


#%% Mean results compare
plt.figure(5)
plt.subplot(221)
plt.plot(Nlevels, NPPaCO2[:,6],'bo-',label='Poly, aCO$_{2}$')
plt.plot(Nlevels, NPPeCO2[:,6],'rD-',label='Poly, eCO$_{2}$')
plt.plot(Nlevels, NPPaCO2[:,3],'bo--',mfc='none',label='$\phi_{RL}$=4, aCO$_{2}$')
plt.plot(Nlevels, NPPeCO2[:,3],'rD--',mfc='none',label='$\phi_{RL}$=4, eCO$_{2}$')
plt.legend()
#plt.xlim((100,450))
plt.ylim((0.1,1.2))
plt.title('(a)', x=0.85, y=0.85,fontdict=font)
#plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('NPP (kgC $m^{-2}$)', fontdict=font)

plt.subplot(222)
plt.plot(Nlevels, NPLaCO2[:,6]/NPPaCO2[:,6],'bo-')
plt.plot(Nlevels, NPLeCO2[:,6]/NPPeCO2[:,6],'rD-')
plt.plot(Nlevels, NPLaCO2[:,3]/NPPaCO2[:,3],'bo--',mfc='none')
plt.plot(Nlevels, NPLeCO2[:,3]/NPPeCO2[:,3],'rD--',mfc='none')
plt.legend()
#plt.xlim((100,450))
plt.ylim((0.08,0.17))
plt.title('(b)', x=0.1, y=0.85,fontdict=font)
#plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
#plt.ylabel('$f_{NPP,leaf}$', fontdict=font)
plt.ylabel('Foliage NPP fraction', fontdict=font)

plt.subplot(223)
plt.plot(Nlevels, NPRaCO2[:,6]/NPPaCO2[:,6],'bo-')
plt.plot(Nlevels, NPReCO2[:,6]/NPPeCO2[:,6],'rD-')
plt.plot(Nlevels, NPRaCO2[:,3]/NPPaCO2[:,3],'bo--',mfc='none')
plt.plot(Nlevels, NPReCO2[:,3]/NPPeCO2[:,3],'rD--',mfc='none')
plt.legend()
#plt.xlim((100,450))
plt.ylim((0,0.75))
plt.title('(c)', x=0.85, y=0.85,fontdict=font)
plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('Root NPP fraction', fontdict=font)

plt.subplot(224)
plt.plot(Nlevels, NPWaCO2[:,6]/NPPaCO2[:,6],'bo-') # ,label='ambient $CO_{2}$'
plt.plot(Nlevels, NPWeCO2[:,6]/NPPeCO2[:,6],'rD-') # ,label='elevated $CO_{2}$'
plt.plot(Nlevels, NPWaCO2[:,3]/NPPaCO2[:,3],'bo--',mfc='none')
plt.plot(Nlevels, NPWeCO2[:,3]/NPPeCO2[:,3],'rD--',mfc='none')
plt.legend()
#plt.xlim((100,450))
plt.ylim((0,0.75))
plt.title('(d)', x=0.1, y=0.85,fontdict=font)
plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('Wood NPP fraction', fontdict=font)

plt.figure(6)
plt.subplot(221)
plt.plot(Nlevels, NPPaCO2[:,6],'bo-',label='Poly, aCO$_{2}$')
plt.plot(Nlevels, NPPeCO2[:,6],'rD-',label='Poly, eCO$_{2}$')
plt.plot(Nlevels, NPPaCO2[:,3],'bo--',mfc='none',label='$\phi_{RL}$=4, aCO$_{2}$')
plt.plot(Nlevels, NPPeCO2[:,3],'rD--',mfc='none',label='$\phi_{RL}$=4, eCO$_{2}$')
plt.legend()
#plt.xlim((100,450))
plt.ylim((0.1,1.2))
plt.title('(a)', x=0.85, y=0.85,fontdict=font)
#plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('NPP (kgC $m^{-2}$)', fontdict=font)

plt.subplot(222)
plt.plot(Nlevels, NPLaCO2[:,6],'bo-')
plt.plot(Nlevels, NPLeCO2[:,6],'rD-')
plt.plot(Nlevels, NPLaCO2[:,3],'bo--',mfc='none')
plt.plot(Nlevels, NPLeCO2[:,3],'rD--',mfc='none')
plt.legend()
#plt.xlim((100,450))
plt.yticks(np.arange(0, 0.15, 0.04))
plt.title('(b)', x=0.1, y=0.85,fontdict=font)
#plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
#plt.ylabel('$f_{NPP,leaf}$', fontdict=font)
plt.ylabel('Foliage NPP', fontdict=font)

plt.subplot(223)
plt.plot(Nlevels, NPRaCO2[:,6],'bo-')
plt.plot(Nlevels, NPReCO2[:,6],'rD-')
plt.plot(Nlevels, NPRaCO2[:,3],'bo--',mfc='none')
plt.plot(Nlevels, NPReCO2[:,3],'rD--',mfc='none')
plt.legend()
#plt.xlim((100,450))
plt.ylim((0,0.4))
plt.title('(c)', x=0.85, y=0.85,fontdict=font)
plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('Root NPP', fontdict=font)

plt.subplot(224)
plt.plot(Nlevels, NPWaCO2[:,6],'bo-') # ,label='ambient $CO_{2}$'
plt.plot(Nlevels, NPWeCO2[:,6],'rD-') # ,label='elevated $CO_{2}$'
plt.plot(Nlevels, NPWaCO2[:,3],'bo--',mfc='none')
plt.plot(Nlevels, NPWeCO2[:,3],'rD--',mfc='none')
plt.legend()
#plt.xlim((100,450))
plt.ylim((0,0.75))
plt.title('(d)', x=0.1, y=0.85,fontdict=font)
plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('Wood NPP', fontdict=font)

#%% Mean results compare
plt.figure(61)
plt.subplot(321)
plt.plot(Nlevels, NPLaCO2[:,6],'bo-',label='Poly, aCO$_{2}$')
plt.plot(Nlevels, NPLeCO2[:,6],'rD-',label='Poly, eCO$_{2}$')
plt.plot(Nlevels, NPLaCO2[:,3],'bo--',mfc='none',label='$\phi_{RL}$=4, aCO$_{2}$')
plt.plot(Nlevels, NPLeCO2[:,3],'rD--',mfc='none',label='$\phi_{RL}$=4, eCO$_{2}$')
plt.legend()
#plt.xlim((100,450))
plt.yticks(np.arange(0, 0.25, 0.05))
plt.title('(a)', x=0.85, y=0.85,fontdict=font)
#plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
#plt.ylabel('$f_{NPP,leaf}$', fontdict=font)
plt.ylabel('Foliage NPP\n(kgC $m^{-2}yr^{-1}$)', fontdict=font)

plt.subplot(323)
plt.plot(Nlevels, NPRaCO2[:,6],'bo-')
plt.plot(Nlevels, NPReCO2[:,6],'rD-')
plt.plot(Nlevels, NPRaCO2[:,3],'bo--',mfc='none')
plt.plot(Nlevels, NPReCO2[:,3],'rD--',mfc='none')
#plt.legend()
#plt.xlim((100,450))
plt.ylim((0,0.8))
plt.title('(c)', x=0.1, y=0.85,fontdict=font)
#plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('Root NPP\n(kgC $m^{-2}yr^{-1}$)', fontdict=font)

plt.subplot(325)
plt.plot(Nlevels, NPWaCO2[:,6],'bo-') # ,label='ambient $CO_{2}$'
plt.plot(Nlevels, NPWeCO2[:,6],'rD-') # ,label='elevated $CO_{2}$'
plt.plot(Nlevels, NPWaCO2[:,3],'bo--',mfc='none')
plt.plot(Nlevels, NPWeCO2[:,3],'rD--',mfc='none')
#plt.legend()
#plt.xlim((100,450))
plt.ylim((0,0.8))
plt.title('(e)', x=0.1, y=0.85,fontdict=font)
plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('Wood NPP\n(kgC $m^{-2}yr^{-1}$)', fontdict=font)

plt.subplot(322)
plt.plot(Nlevels, NPLaCO2[:,6]/NPPaCO2[:,6],'bo-')
plt.plot(Nlevels, NPLeCO2[:,6]/NPPeCO2[:,6],'rD-')
plt.plot(Nlevels, NPLaCO2[:,3]/NPPaCO2[:,3],'bo--',mfc='none')
plt.plot(Nlevels, NPLeCO2[:,3]/NPPeCO2[:,3],'rD--',mfc='none')
#plt.legend()
#plt.xlim((100,450))
plt.yticks(np.arange(0, 0.25, 0.05))
plt.title('(b)', x=0.1, y=0.85,fontdict=font)
#plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
#plt.ylabel('$f_{NPP,leaf}$', fontdict=font)
plt.ylabel('Foliage NPP fraction', fontdict=font)

plt.subplot(324)
plt.plot(Nlevels, NPRaCO2[:,6]/NPPaCO2[:,6],'bo-')
plt.plot(Nlevels, NPReCO2[:,6]/NPPeCO2[:,6],'rD-')
plt.plot(Nlevels, NPRaCO2[:,3]/NPPaCO2[:,3],'bo--',mfc='none')
plt.plot(Nlevels, NPReCO2[:,3]/NPPeCO2[:,3],'rD--',mfc='none')
#plt.legend()
#plt.xlim((100,450))
plt.ylim((0,0.8))
plt.title('(d)', x=0.1, y=0.85,fontdict=font)
#plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('Root NPP fraction', fontdict=font)

plt.subplot(326)
plt.plot(Nlevels, NPWaCO2[:,6]/NPPaCO2[:,6],'bo-') # ,label='ambient $CO_{2}$'
plt.plot(Nlevels, NPWeCO2[:,6]/NPPeCO2[:,6],'rD-') # ,label='elevated $CO_{2}$'
plt.plot(Nlevels, NPWaCO2[:,3]/NPPaCO2[:,3],'bo--',mfc='none')
plt.plot(Nlevels, NPWeCO2[:,3]/NPPeCO2[:,3],'rD--',mfc='none')
#plt.legend()
#plt.xlim((100,450))
plt.ylim((0,0.8))
plt.title('(f)', x=0.1, y=0.85,fontdict=font)
plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('Wood NPP fraction', fontdict=font)


#%%
plt.figure(7)
plt.subplot(221)
plt.plot(Nlevels, BMaCO2[:,6],'bo-',label='Poly, aCO$_{2}$')
plt.plot(Nlevels, BMeCO2[:,6],'rD-',label='Poly, eCO$_{2}$')
plt.plot(Nlevels, BMaCO2[:,4],'bo--',mfc='none',label='$\phi_{RL}$=4, aCO$_{2}$')
plt.plot(Nlevels, BMeCO2[:,4],'rD--',mfc='none',label='$\phi_{RL}$=4, eCO$_{2}$')
plt.legend()
plt.yticks(np.arange(0, 55, 10))
plt.title('(a)', x=0.85, y=0.85,fontdict=font)
#plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('Plant biomass (kgC m$^{-2}$)', fontdict=font)

plt.subplot(222)
plt.plot(Nlevels, BMeCO2[:,6]-BMaCO2[:,6],'bo-',label='Poly')
plt.plot(Nlevels, BMeCO2[:,4]-BMaCO2[:,3],'bo--',mfc='none',label='$\phi_{RL}$=4')
plt.legend()
plt.yticks(np.arange(-4, 20, 4))
plt.title('(b)', x=0.7, y=0.85,fontdict=font)
#plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('Difference (kgC m$^{-2}$)', fontdict=font)

plt.subplot(223)
plt.plot(Nlevels, BMaCO2[:,0:6],'o--',mfc='none')
plt.plot(Nlevels, BMaCO2[:,6],'bo-')
#plt.xlim((100,450))
plt.ylim((0,50))
plt.title('(c)', x=0.8, y=0.85,fontdict=font)
plt.xlabel('Total N (gN m$^{-2}$)', fontdict=font)
plt.ylabel('Plant biomass (kgC $m^{-2}$)', fontdict=font)
plt.legend( ('$\phi_{RL}$: 1', '$\phi_{RL}$: 2', '$\phi_{RL}$: 3', \
             '$\phi_{RL}$: 4', '$\phi_{RL}$: 5', '$\phi_{RL}$: 6') \
            , ncol=2, loc=2) # , bbox_to_anchor=(1, 0.5) 

plt.subplot(224)
plt.plot(Nlevels, BMeCO2[:,0:6],'o--',mfc='none')
plt.plot(Nlevels, BMeCO2[:,6],'rD-')
#plt.xlim((100,450))
plt.ylim((0,50))

plt.title('(d)', x=0.15, y=0.85,fontdict=font)
#plt.title('Plant Biomass', x=0.8, y=0.85,fontdict=font)
plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('Plant biomass (kgC $m^{-2}$)', fontdict=font)
#%%
#Critical height
plt.figure(8)
plt.plot(Nlevels, HTaCO2[:,6],'bo-',label='Poly, aCO$_{2}$')
plt.plot(Nlevels, HTeCO2[:,6],'rD-',label='Poly, eCO$_{2}$')
plt.plot(Nlevels, HTaCO2[:,3],'bo--',mfc='none',label='$\phi_{RL}$=4, aCO$_{2}$')
plt.plot(Nlevels, HTeCO2[:,3],'rD--',mfc='none',label='$\phi_{RL}$=4, eCO$_{2}$')
plt.legend()
#plt.title('(b)', x=0.85, y=0.85,fontdict=font)
plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.ylabel('Critical Height (m)', fontdict=font)

#%% Plot LAI
plt.figure(10)
plt.subplot(211)
plt.plot(Nlevels,(LAIaCO2[:,0:6]),'o--',mfc='none')
plt.plot(Nlevels,(LAIaCO2[:,6]),'bo-')
#plt.title('(f)', x=0.15, y=0.85,fontdict=font)
plt.ylim((0.,4))
plt.legend( ('$\phi_{RL}$: 1', '$\phi_{RL}$: 2', '$\phi_{RL}$: 3', \
             '$\phi_{RL}$: 4', '$\phi_{RL}$: 5', '$\phi_{RL}$: 6', 'poly'),\
              loc='best',ncol=2)
plt.title('(a): aCO$_{2}$', x=0.8, y=0.85,fontdict=font)
plt.ylabel('Leaf area index', fontdict=font)
#plt.xlabel('Total N ($gN/m^{2}$)', fontdict=font)

plt.subplot(212)
plt.plot(Nlevels,(LAIeCO2[:,0:6]),'o--',mfc='none')
plt.plot(Nlevels,(LAIeCO2[:,6]),'rD-')
#plt.title('(f)', x=0.15, y=0.85,fontdict=font)
plt.ylim((0.,4))
plt.legend( ('$\phi_{RL}$: 1', '$\phi_{RL}$: 2', '$\phi_{RL}$: 3', \
             '$\phi_{RL}$: 4', '$\phi_{RL}$: 5', '$\phi_{RL}$: 6', 'poly'),\
              loc='best',ncol=2)
plt.title('(b): eCO$_{2}$', x=0.8, y=0.85,fontdict=font)
plt.ylabel('Leaf area index', fontdict=font)
plt.xlabel('Total N ($gN/m^{2}$)', fontdict=font)

plt.figure(101)
plt.subplot(221)
plt.plot(Nlevels,(LAIaCO2[:,0:6]),'o--',mfc='none')
plt.plot(Nlevels,(LAIaCO2[:,6]),'bo-')
#plt.title('(f)', x=0.15, y=0.85,fontdict=font)
plt.ylim((0.,4.5))
plt.legend( ('$\phi_{RL}$: 1', '$\phi_{RL}$: 2', '$\phi_{RL}$: 3', \
             '$\phi_{RL}$: 4', '$\phi_{RL}$: 5', '$\phi_{RL}$: 6'),\
              loc='best',ncol=2)
plt.title('(a): aCO$_{2}$', x=0.8, y=0.85,fontdict=font)
plt.ylabel('Leaf area index', fontdict=font)
#plt.xlabel('Total N ($gN/m^{2}$)', fontdict=font)

plt.subplot(222)
plt.plot(Nlevels,(LAIeCO2[:,0:6]),'o--',mfc='none')
plt.plot(Nlevels,(LAIeCO2[:,6]),'rd-')
#plt.title('(f)', x=0.15, y=0.85,fontdict=font)
plt.ylim((0.,4.5))
plt.legend( ('$\phi_{RL}$: 1', '$\phi_{RL}$: 2', '$\phi_{RL}$: 3', \
             '$\phi_{RL}$: 4', '$\phi_{RL}$: 5', '$\phi_{RL}$: 6'),\
              loc=4,ncol=2)
plt.title('(b): eCO$_{2}$', x=0.8, y=0.85,fontdict=font)
plt.ylabel('Leaf area index', fontdict=font)
plt.xlabel('Total N ($gN/m^{2}$)', fontdict=font)

plt.subplot(223)
plt.plot(Nlevels, 100*(LAIeCO2[:,3]-LAIaCO2[:,3])/LAIaCO2[:,3],'ro-',)
plt.legend('Mono',loc='best')
#plt.yticks(np.arange(-60, 0, 10))
plt.title('(c)', x=0.8, y=0.85,fontdict=font)
plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.tick_params(labelbottom='on')
plt.ylabel('Relative changes (%)', fontdict=font)

plt.subplot(224)
plt.plot(Nlevels, 100*(LAIeCO2[:,6]-LAIaCO2[:,6])/LAIaCO2[:,6],'rs-')
plt.legend('Poly',loc=4)
#plt.yticks(np.arange(-60, 0, 10))
plt.title('(d)', x=0.8, y=0.85,fontdict=font)
plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)
plt.tick_params(labelbottom='on')
plt.ylabel('Relative changes (%)', fontdict=font)

#%% Nitrogen patterns
plt.figure(11)
plt.subplot(211)
plt.plot(Nlevels[0:6],NmaCO2[0:6,0:6],'o',linestyle = 'None')
plt.title('(a)', x=0.85, y=0.9,fontdict=font)
#plt.xlim((0,6))
plt.ylim((0.05,0.2))
plt.ylabel('mineral N\n(gN $m^{-2}$)', fontdict=font)
plt.xlabel('Total N ($gN/m^{2}$)', fontdict=font)
plt.legend( ('$\phi_{RL}$: 1', '$\phi_{RL}$: 2', '$\phi_{RL}$: 3', \
             '$\phi_{RL}$: 4', '$\phi_{RL}$: 5', '$\phi_{RL}$: 6', \
             '$\phi_{RL}$: 7', '$\phi_{RL}$: 8'), \
             loc='best', ncol=2)

plt.subplot(212)
plt.plot(Nlevels[0:6],NmeCO2[0:6,0:6],'o',linestyle = 'None')
plt.title('(b)', x=0.85, y=0.9,fontdict=font)
plt.ylim((0.05,0.2))
plt.ylabel('mineral N\n(gN $m^{-2}$)', fontdict=font)
plt.xlabel('Total N (gN $m^{-2}$)', fontdict=font)

plt.figure(12)
plt.plot(Nlevels[0:6], NmaCO2[0:6,6],'bo-',label='aCO$_{2}$')
plt.plot(Nlevels[0:6], NmeCO2[0:6,6],'rD-',label='eCO$_{2}$')
plt.legend()
plt.ylim((0.09,0.15))
#plt.title('(b)', x=0.85, y=0.85,fontdict=font)
plt.xlabel('Total N (gN m$^{-2}$)', fontdict=font)
plt.ylabel('Mineral N (gN m$^{-2}$)', fontdict=font)

plt.figure(13)
plt.subplot(211)
plt.plot(Nlevels[0:6], NmaCO2[0:6:,0:6],'o')
plt.plot(Nlevels[0:6], NmaCO2[0:6:,6],'bD-',label='aCO$_{2}$')
#plt.xlim((100,450))
#plt.ylim((0.15,0.35))
plt.ylim((0.08,0.2))
plt.legend( ('$\phi_{RL}$: 1', '$\phi_{RL}$: 2', '$\phi_{RL}$: 3', \
             '$\phi_{RL}$: 4', '$\phi_{RL}$: 5', '$\phi_{RL}$: 6', 'Poly') \
            , ncol=1, loc='center left', bbox_to_anchor=(1, 0.5) )
plt.title('(a): aCO$_{2}$', x=0.15, y=0.85,fontdict=font)
#plt.xlabel('Total N (gN m$^{-2}$)', fontdict=font)
plt.ylabel('Mineral N (gN m$^{-2}$)', fontdict=font)

plt.subplot(212)
plt.plot(Nlevels[0:6], NmeCO2[0:6,0:6],'o')
plt.plot(Nlevels[0:6], NmeCO2[0:6,6],'bD-',label='eCO$_{2}$')
#plt.xlim((100,450))
plt.ylim((0.08,0.2))
#plt.legend( ('$\phi_{RL}$: 1', '$\phi_{RL}$: 2', '$\phi_{RL}$: 3', \
#             '$\phi_{RL}$: 4', '$\phi_{RL}$: 5', '$\phi_{RL}$: 6','Multi') )
plt.title('(b): eCO$_{2}$', x=0.15, y=0.85,fontdict=font)
plt.xlabel('Total N (gN m$^{-2}$)', fontdict=font)
plt.ylabel('Mineral N (gN m$^{-2}$)', fontdict=font)
