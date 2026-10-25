#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 07 21:52, 2025

Workstation IP: 192.168.0.31

@author: eweng
"""
import gc
import numpy as np
import netCDF4 as nc4
import  os
import sys
from datetime import datetime
import gzip
import csv
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler

linestyle_cycler = cycler('linestyle',['-','-.','--',':'])

font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 14,
        }
PI = 3.1415926
# CRU data parameters
N0_Lat = 360 # Dimension of CRU data (0.5x0.5)
N0_Lon = 720

PFT_ID = ['C4G',  'C3G',  'TrE',  'TrD',  'TmE',  'TmD',  'Nfx',  'DeS']
PFTnames = ['C4 grass',\
            'C3 grass',\
            'Tropical evergreen',\
            'Tropical deciduous',\
            'Temperate evergreen',\
            'Cold deciduous',\
            'Nitrogen fixers',\
            'Desert shrubs']

PFTUnit = ['m2/m2','m2/m2','m2/m2','m2/m2','m2/m2','m2/m2','m2/m2','m2/m2']

# Ecosystem file variables
EcoVars = ['CAI',    'LAI',  'GPP', 'Rauto', 'Rh', 'Burned', \
           'Tavg', 'Rain', 'SoilWater','Transp','Evap','Runoff',\
           'plantC', 'soilC', 'plantN','soilN', \
           'NSC',    'SeedC','leafC', 'rootC', 'swC', 'hwC',\
           'NSN',    'SeedN','leafN', 'rootN', 'swN', 'hwN',\
           'fineL',  'strucL', 'McrbC', 'fastSOC', 'slowSOC',\
           'fineN',  'strucN', 'McrbN', 'fastSON', 'slowSON', 'mineralN',\
           'WC1_5',  'WC2_25', 'WC3_50', 'WC4_100', 'WC5_120', \
           'N_fxed', 'N_uptk', 'Nm_SL', 'Nm_FR','dNorg','dNgas','dNmin', \
           'treeCA', 'grassCA', 'BMgrass','PET', 'Frisk', 'Pburn','CH4',\
           'HarvC','HarvN','mu','muC','Indv']

EcoLongID = ['Crown area index','Leaf area index','Gross Primary Production',\
             'Autotrophic respiration','Heterotrophic Respiration','Burned carbon',\
            'Yearly mean temperature','Yearly rainfall','Soil water amount',\
            'Yearly transpiration','Yearly evaporation','Yearly runoff',\
            'Plant Biomass','Soil Organic Matter','Plant nitrogen','Soil nitrogen',\
            'Non-structural carbon','Seed carbon', 'leaf carbon', 'root carbon', \
            'Sapwood carbon', 'Heartwood carbon',\
            'Non-structural nitrogen','Seed nitrogen', 'leaf nitrogen', 'root nitrogen', \
            'Sapwood nitrogen', 'Heartwood nitrogen',\
            'fineL',  'strucL', 'McrbC', 'fastSOC', 'slowSOC',\
            'fineN',  'strucN', 'McrbN', 'fastSON', 'slowSON', 'mineralN',\
            'WC1_5',  'WC2_25', 'WC3_50', 'WC4_100', 'WC5_120', \
            'N_fixed', 'N uptake', 'Nm_SL', 'Nm_FR','Nloss1','Nloss2','Nloss3',\
            'Woody crown area','Grass crown area', 'BMgrass', \
            'Potential evapotranspiration','Fire risk','Fire probability', 'Methane', \
            'Harvested carbon','Harvested nitrogen',\
            'Mortality rate','Mortality carbon flux','Woody individuals']

EcoUnit = ['m2/m2','m2/m2','KgC m-2 yr-1','KgC m-2 yr-1','KgC m-2 yr-1','KgC m-2 yr-1',\
           'degree C','mm/year','mm','mm/year','mm/year','mm/year',\
           'KgC m-2','KgC m-2','gN m-2','gN m-2',\
           'KgC m-2','KgC m-2','KgC m-2','KgC m-2','KgC m-2','KgC m-2',\
           'gN m-2','gN m-2','gN m-2','gN m-2','gN m-2','gN m-2',\
           'KgC m-2','KgC m-2','KgC m-2','KgC m-2','KgC m-2',\
           'gN m-2','gN m-2','gN m-2','gN m-2','gN m-2','gN m-2',\
           'mm','mm','mm','mm','mm','gN m-2 yr-1','gN m-2 yr-1',\
           'gN m-2 yr-1','gN m-2 yr-1','gN m-2 yr-1','gN m-2 yr-1','gN m-2 yr-1',\
           'm2/m2','m2/m2','KgC m-2','mm/year','times/yr','times/yr','KgC m-2 yr-1', \
           'KgC m-2 yr-1','gN m-2 yr-1','fraction yr-1','KgC m-2 yr-1','individuals/m2']

subfolder = ['Ecosystem','Cohort']

N_pfts   = 8 # total PFTs at one site, 4
N_Layers = 3
Npre     = 6
N_gridV  = 0 # 58 # added two more N losses and CH4, 55 # 54, added YealryTmp, 11/23/2025
totYrs   = 0 # Will be updated by reading an ecosystem data file


def read_yearly_csv(fname_gz, ncols):
    """Read the first ncols columns of a gzipped yearly csv (header skipped) as float64.
    Much lighter than np.array(list(csv.reader(f))), which holds every value as a Python
    string and exhausted memory over tens of thousands of files."""
    return pd.read_csv(fname_gz, compression='gzip', header=0, skipinitialspace=True,
                       usecols=range(ncols), dtype=np.float64).to_numpy()

#%% --------------------- Check the files ---------------------
# 'eCO2' # 'N2g1123' #  'N3g1121' # 'N3gLowNfx' # 'N3gTr10' # 'Ndps3g'
# 'N4g1128' #  'Warming2C' # '0.5LonLat_N2g1125'
expID = 'N2gCO21901' # 'N2gLowNout' # 'BaseN2gThnG' # 'N3gWmu0Low' # 'TmIgnN3g' # 'MI0Fr2N3g' 
#Resolution = 1 #Grids selected
# /media/eweng/HD2/weng/GlobalESSPFTs/Simulations/N2gCO21901/Run_N2gCO21901
path0 = '/media/eweng/HD2/weng/GlobalESSPFTs/Simulations/GlobalRun_'+expID+'/'
path0 = '/media/eweng/HD2/weng/GlobalESSPFTs/Simulations/' + expID + '/Run_' + expID+ '/'
fpath = path0 +'Yearly/'
fpout = path0 +''

ecofiles = [f for f in os.listdir(fpath) if "Ecosystem_yearly.csv.gz" in f]
cohfiles = [f for f in os.listdir(fpath) if "Cohort_yearly.csv.gz" in f]

ecofiles.sort()
cohfiles.sort()

N_files = len(ecofiles)
if len(ecofiles) != len(cohfiles):
    print(f"Ecosystem and Cohort files are not equal: {len(ecofiles)} != {len(cohfiles)}")
    sys.exit("Script terminated because values are not equal.")
else:
    print(f"Ecosystem and Cohort files are equal: {len(ecofiles)} == {len(cohfiles)}")


#%% --------------------- Find out Resolution, and file data lines and rows
GridID    = np.zeros((N_files), dtype='int')
LonFiles  = np.zeros((N_files), dtype='int')
LatFiles  = np.zeros((N_files), dtype='int')
# Find out iLon and iLat
maxLat = -1
minLat = 999
maxLon = -1
minLon = 999
for ifile in range(N_files):
    # Read the file names to get Lat and Lon
    grid_str = ecofiles[ifile][Npre:Npre+6]
    iLonLat = int(grid_str)
    GridID[ifile]  = iLonLat
    #print (iLonLat)
    iLon = int(iLonLat/1000)
    iLat  = iLonLat - iLon*1000
    LonFiles[ifile]  = iLon
    LatFiles[ifile]  = iLat
    maxLat = max(maxLat,iLat)
    minLat = min(minLat,iLat)
    maxLon = max(maxLon,iLon)
    minLon = min(minLon,iLon)

# Get model run years and the number of ecosystem variables
fname_gz = fpath + ecofiles[1]
try:
    with gzip.open(fname_gz, 'rt', newline='') as f:
        reader = csv.reader(f)
        header = next(reader)  # Read the header row
        #print(f"Header: {header}")
        LandYrV=(np.array(list(csv.reader(f,delimiter=','))))
except IOError as e:
    print(f"An I/O error occurred: {e}")

f.close()

totYrs  = len(LandYrV)
N_gridV = len(LandYrV[1]) - 1 - 2 # two empty and one additional cols in each row
YR0 = int(totYrs * 0.2) # for calculating averages

# Obtain the resolution of model run
Resolution = int(GridID[1] -GridID[0])
r_list = [1, 2, 3, 4, 5]
if Resolution in r_list:
    print(f"Grid resolution is {Resolution}.")
    print(N_files,N_gridV, totYrs)
    print(minLat,maxLat,minLon,maxLon)
else:
    print(f"Resolution {Resolution} is not in the resolution list.")
    sys.exit(0)

#N_Lat = int((maxLat - minLat)/Resolution) + 1
#N_Lon = int((maxLon - minLon)/Resolution) + 1
N_Lat = int(N0_Lat/Resolution)
N_Lon = int(N0_Lon/Resolution)

#%% --------------------- Check simulated grids coverage ---------------------
SimuCoverage = np.zeros((N_Lat,N_Lon))
for ifile in range(N_files):
    iLon = int(LonFiles[ifile]/Resolution)-1
    iLat = int(LatFiles[ifile]/Resolution)-1
    SimuCoverage[iLat,iLon] = 1.0

fig = plt.figure(1)
plt.clf()
plt.title('Simulated Grids')
plt.imshow(np.flipud(SimuCoverage),interpolation='none') #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)


#%% --------------------- Read in ecosystem files ---------------------
# Opent output files
feco  = open(fpout + "EcoFileNames.txt", "w")

# For ecosystem data
AvgGridsData = np.zeros((N_gridV+3, N_Lat, N_Lon)) # include mu and muC
VldFiles     = np.zeros((N_files), dtype='int')
N_valid      = 0
for ifile in range(N_files):
    iLonLat = GridID[ifile]
    iLon = LonFiles[ifile] - 1
    iLat = LatFiles[ifile] - 1
    m = int(iLon/Resolution)
    n = int(iLat/Resolution)

    # Open and read ecosystem file
    fname_gz = fpath + ecofiles[ifile]
    try:
        LandYr = read_yearly_csv(fname_gz, 2 + N_gridV)
        VldFiles[ifile] = 1
        N_valid = N_valid + 1
    except (IOError, EOFError, ValueError) as e:
        print(f"Cannot read {fname_gz}: {e}")
        continue # skip this cycle

    feco.write(ecofiles[ifile] + '\n')

    # Skip files with too few years
    rows = LandYr.shape[0]
    if rows < totYrs/10:
        continue
    if N_valid % 1000 == 0: print (N_valid,ecofiles[ifile])
    AvgGridsData[0:N_gridV, n, m] = np.mean(LandYr[YR0:rows,2:2+N_gridV], axis=0)
    del LandYr

feco.close()

#%% --------------------- Read in cohort files ---------------------
# Opent output files
fcoh  = open(fpout + "CohFileNames.txt", "w")
fcai  = open(fpout + "CA_Tmp.txt", 'w')
#fden  = open(fpout + "Den_Tmp.txt", 'w')
#fmu   = open(fpout + "Mu_Tmp.txt", 'w')

# For cohort data
meanPFTGPP = np.zeros((N_pfts,N_Lat,N_Lon))
meanPFTNPP = np.zeros((N_pfts,N_Lat,N_Lon))
meanPFTBA  = np.zeros((N_pfts,N_Lat,N_Lon))
meanPFTCA  = np.zeros((N_pfts,N_Lat,N_Lon))
meanPFTLA  = np.zeros((N_pfts,N_Lat,N_Lon))
meanPFTBM  = np.zeros((N_pfts,N_Lat,N_Lon))
meanPFTHT  = np.zeros((N_pfts,N_Lat,N_Lon))
meanPFTden = np.zeros((N_pfts,N_Lat,N_Lon))
meanPFTmu  = np.zeros((N_pfts,N_Lat,N_Lon))
meanPFTmuC = np.zeros((N_pfts,N_Lat,N_Lon))

GPP = np.zeros((totYrs, N_pfts))
NPP = np.zeros((totYrs, N_pfts))
BA  = np.zeros((totYrs, N_pfts))
CA  = np.zeros((totYrs, N_pfts))
LA  = np.zeros((totYrs, N_pfts))
BM  = np.zeros((totYrs, N_pfts))
HT  = np.zeros((totYrs, N_pfts))
den = np.zeros((totYrs, N_pfts))
mu  = np.zeros((totYrs, N_pfts))
muC = np.zeros((totYrs, N_pfts))

# Read in cohort files
N_valid = 0
for ifile in range(N_files):
    iLonLat = GridID[ifile]
    iLon = LonFiles[ifile] - 1
    iLat = LatFiles[ifile] - 1
    m = int(iLon/Resolution)
    n = int(iLat/Resolution)

    # Open and read Cohort file
    fname_gz = fpath + cohfiles[ifile]
    try:
        CCYr = read_yearly_csv(fname_gz, 41)
        N_valid = N_valid + 1
    except (IOError, EOFError, ValueError) as e:
        print(f"Cannot read {fname_gz}: {e}")
        continue # skip this cycle
    if N_valid % 1000 == 0: print (N_valid,cohfiles[ifile])

    # Calculate PFT-level  GPP, NPP, BA, CA, BM, LAI, height (vectorized over cohorts)
    iYr    = CCYr[:,1].astype(np.int64) - 1
    iPFT   = CCYr[:,4].astype(np.int64)
    top    = np.minimum(2, CCYr[:,5].astype(np.int64) - 1) == 0 # Only the top layer
    keep   = (iYr >= 0) & (iYr < totYrs) & (iPFT >= 0) & (iPFT < N_pfts)
    idx    = iYr*N_pfts + iPFT
    nTot   = totYrs*N_pfts
    Ndens  = CCYr[:,6]
    bmass  = np.sum(CCYr[:,15:21], axis=1)

    def PFTsum(w, sel):
        ok = keep & sel
        return np.bincount(idx[ok], weights=w[ok], minlength=nTot).reshape(totYrs, N_pfts)

    allc = np.ones(len(CCYr), dtype=bool)
    GPP = PFTsum(Ndens*CCYr[:,22]/10000, allc)
    NPP = PFTsum(Ndens*CCYr[:,23]/10000, allc)
    BA  = PFTsum(Ndens*PI*0.25*CCYr[:,11]**2, allc)
    LA  = PFTsum(Ndens*CCYr[:,14]/10000, allc)
    BM  = PFTsum(Ndens*bmass/10000, allc)
    HT  = np.zeros(nTot)
    np.maximum.at(HT, idx[keep], CCYr[keep,12])
    HT  = HT.reshape(totYrs, N_pfts)
    den = PFTsum(Ndens, top)
    muW = PFTsum(Ndens*CCYr[:,29], top)
    muC = PFTsum(Ndens*CCYr[:,29]*bmass/10000, top)
    CA  = PFTsum(Ndens*CCYr[:,13]/10000, top)
    # Density weighted mortality rate
    mu  = np.divide(muW, den, out=np.zeros_like(muW), where=den > 1e-4)

    # Mean
    #iLat 175, iLong 251
    meanPFTGPP[:,n,m] = np.mean(GPP[YR0:totYrs,:],axis=0)
    meanPFTNPP[:,n,m] = np.mean(NPP[YR0:totYrs,:],axis=0)
    meanPFTmuC[:,n,m] = np.mean(muC[YR0:totYrs,:],axis=0)
    meanPFTBA [:,n,m] = np.mean(BA [YR0:totYrs,:],axis=0)
    meanPFTCA [:,n,m] = np.mean(CA [YR0:totYrs,:],axis=0)
    meanPFTLA [:,n,m] = np.mean(LA [YR0:totYrs,:],axis=0)
    meanPFTBM [:,n,m] = np.mean(BM [YR0:totYrs,:],axis=0)
    meanPFTHT [:,n,m] = np.mean(HT [YR0:totYrs,:],axis=0)
    meanPFTden[:,n,m] = np.mean(den[YR0:totYrs,:],axis=0)
    for k in range(8):
        if meanPFTden[k,n,m] > 1.0:
            totMu = 0.0
            totDn = 0.0
            for j in range(YR0,totYrs):
                totMu = totMu + den[j, k] * mu[j, k]
                totDn = totDn + den[j, k]
            meanPFTmu [k, n, m] = totMu/totDn
    
    #print('Mortality rate:')
    #print(meanPFTmu [:, n, m])


    # Remove CCYr variables and release memory
    del CCYr

    # --------------- Write out Grid temporal files ---------------
    fcoh.write(cohfiles[ifile] + '\n')
    # Write temporal CA of the PFTs
    for i in range(N_pfts):
        formatted_row = [f"{num:.2f}" for num in CA[:,i]]
        fcai.write(",".join(formatted_row) + "\n")
'''
    # Write temporal mu of the PFTs
    for i in range(N_pfts):
        formatted_row = [f"{num:.2e}" for num in mu[:,i]]
        fmu.write(",".join(formatted_row) + "\n")

    # Write temporal density of the PFTs
    for i in range(N_pfts):
        formatted_row = [f"{num:.2f}" for num in den[:,i]]
        fden.write(",".join(formatted_row) + "\n")
fmu.close()
fden.close()
'''
fcoh.close()
fcai.close()


#%% ----------------- Write data to files --------------------------

#Calculate mu and muC, and write txt output for all grids
# Output files
fca   = open(fpout + "Eco_CA_PFTs.txt", 'w')
fgrid = open(fpout + "GridIDs.txt", "w")

# Calculate Woody PFT's density, mu and muC
AvgMu = np.zeros((N_Lat, N_Lon))
tmp   = np.zeros((totYrs, N_pfts))

SumMuC    = np.nansum(meanPFTmuC, axis=0)
woody_den = meanPFTden[2:7, :, :]
woody_mu  = meanPFTmu [2:7, :, :]
SumDen    = np.nansum(woody_den, axis=0) # Trees only

# Density-weighted mean of meanPFTmu for woody PFTs
valid = (woody_den > 0.0) & (woody_mu > 0.0)
weighted_mu = np.where(valid, woody_den * woody_mu, 0.0)
np.divide(np.sum(weighted_mu, axis=0), SumDen, out=AvgMu, where=SumDen > 1.0e-6)
for ifile in range(N_files):
    iLonLat = GridID[ifile]
    m = int((LonFiles[ifile] - 1) / Resolution)
    n = int((LatFiles[ifile] - 1) / Resolution)

    fgrid.write(f"{iLonLat}" + '\n')

    # Write men CA and GPP to csv file
    # Tmp,Rain,PET,LAI,GPP,BM,Frisk,Pburn, PFTs (8 CA)
    #formatted_row = [f"{num:.2f}" for num in meanPFTCA [:,n,m]]
    #row = np.concatenate((AvgGridsData[[6,51,1,2,11,52,53],n,m],
    row = np.concatenate((AvgGridsData[[6,7,54,1,2,12,55,56],n,m],
          meanPFTCA[:, n, m], meanPFTmu [:,n,m], meanPFTden[:,n,m]))
    formatted_row = [f"{x:.4f}" for x in row]
    fca.write(",".join(formatted_row) + "\n")

fca.close()
fgrid.close()

#% Write to netCDF files
LonLatStep = 0.5 * Resolution
Lon0  = -180.0 + LonLatStep * 0.5
Lat0  = -90.0  + LonLatStep * 0.5
today = datetime.now()
Nvars = len(EcoVars)

# Create a netcdf dataset for ecosystem variables
f1 = nc4.Dataset(fpout + 'BiomeE_Simu_' + expID + '.nc','w', format='NETCDF4') #'w' write
f1.createDimension('lon', N_Lon)
f1.createDimension('lat', N_Lat)
longitude = f1.createVariable('Longitude', 'f4', 'lon')
latitude  = f1.createVariable('Latitude',  'f4', 'lat')

longitude[:] = np.arange(Lon0,180.0,LonLatStep) #lon
latitude[:]  = np.arange(Lat0,90.0, LonLatStep) #lat

# Put AvgMu to AvgGridsData
AvgGridsData[N_gridV,:,:]   = AvgMu[:,:]
AvgGridsData[N_gridV+1,:,:] = SumMuC[:,:]
AvgGridsData[N_gridV+2,:,:] = SumDen[:,:]

# Write out ecosystem variables
NCout=[]
for iVar in range(Nvars):
        NCout.append(f1.createVariable(EcoVars [iVar], 'f4', ('lat','lon'), zlib=True, complevel=4, shuffle=True))
        NCout[iVar][:,:]      = AvgGridsData[iVar,:,:]
        NCout[iVar].long_name = EcoLongID[iVar]
        NCout[iVar].units     = EcoUnit[iVar]

# Write out PFT distribution
PFTout=[]
for iPFT in range(N_pfts):
        PFTout.append(f1.createVariable(PFT_ID[iPFT], 'f4', ('lat','lon'), zlib=True, complevel=4, shuffle=True))
        PFTout[iPFT].long_name = PFTnames[iPFT]
        PFTout[iPFT].units     = PFTUnit[iPFT]
        PFTout[iPFT][:,:]      = meanPFTCA[iPFT,:,:]

# Write out PFT mortality
PFTu_ID = ['uC4G',  'uC3G',  'uTrE',  'uTrD',  'uTmE',  'uTmD',  'uNfx',  'uDeS']
MuUnit = ['yr-1','yr-1','yr-1','yr-1','yr-1','yr-1','yr-1','yr-1']
Mort=[]
for iPFT in range(N_pfts):
        Mort.append(f1.createVariable(PFTu_ID[iPFT], 'f4', ('lat','lon'), zlib=True, complevel=4, shuffle=True))
        Mort[iPFT].long_name = PFTnames[iPFT]
        Mort[iPFT].units     = MuUnit[iPFT]
        Mort[iPFT][:,:]      = meanPFTmu[iPFT,:,:]

# Write out PFT individual density
PFTd_ID = ['dC4G',  'dC3G',  'dTrE',  'dTrD',  'dTmE',  'dTmD',  'dNfx',  'dDeS']
DenUnit = ['ha-1','ha-1','ha-1','ha-1','ha-1','ha-1','ha-1','ha-1']
Qpft=[]
for iPFT in range(N_pfts):
        Qpft.append(f1.createVariable(PFTd_ID[iPFT], 'f4', ('lat','lon'), zlib=True, complevel=4, shuffle=True))
        Qpft[iPFT].long_name = PFTnames[iPFT]
        Qpft[iPFT].units     = DenUnit[iPFT]
        Qpft[iPFT][:,:]      = meanPFTden[iPFT,:,:]
        
f1.description = "BiomeE simulated variables in the ecosystem output files and PFT distribution"
f1.history = "Created " + today.strftime("%d/%m/%y")

f1.close()

#% Create a netcdf dataset for PFT crown area fraction
N_pfts = len(PFT_ID)
f1 = nc4.Dataset(fpout + 'BiomeE_simulated_PFT_distribution' + '.nc','w', format='NETCDF4') #'w' write
f1.createDimension('lon', N_Lon)
f1.createDimension('lat', N_Lat)
#f1.createDimension('pft', N_pfts)

longitude = f1.createVariable('Longitude', 'f4', 'lon')
latitude  = f1.createVariable('Latitude',  'f4', 'lat')

longitude[:] = np.arange(Lon0,180.0,LonLatStep) #lon
latitude[:]  = np.arange(Lat0,90.0, LonLatStep) #lat

PFTout=[]
for iPFT in range(N_pfts):
        PFTout.append(f1.createVariable(PFT_ID[iPFT], 'f4', ('lat','lon'), zlib=True, complevel=4, shuffle=True))
        PFTout[iPFT].long_name = PFTnames[iPFT]
        PFTout[iPFT].units     = PFTUnit[iPFT]
        PFTout[iPFT][:,:]      = meanPFTCA[iPFT,:,:]
        
f1.description = "PFTs: \
      0:C4G,  1:C3G,  2:TrE,  3:TrD,  4:TmE,  5:TmD,  6:Nfx,  7:DeS "
f1.history = "Created " + today.strftime("%d/%m/%y")

f1.close()


#%% --------------------- Data analysis ---------------------
# Copy to specific variables
meanRain       = AvgGridsData[7,:,:]
meanGPP        = AvgGridsData[2,:,:]  # np.mean(AllGridsData[2,200:totYrs,:,:])
meanRsp        = AvgGridsData[3,:,:]
meanPlantC     = AvgGridsData[12,:,:] # np.mean(AllGridsData[11,200:totYrs,:,:])
meanSoilC      = AvgGridsData[13,:,:]
meanPlantN     = AvgGridsData[14,:,:] # np.mean(AllGridsData[11,200:totYrs,:,:])
meanSoilN      = AvgGridsData[15,:,:]
meanCAtree     = AvgGridsData[51,:,:] # np.mean(AllGridsData[48,200:totYrs,:,:])
meanCAgrass    = AvgGridsData[52,:,:] # np.mean(AllGridsData[49,200:totYrs,:,:])
meanFrisk      = AvgGridsData[55,:,:] # np.mean(AllGridsData[52,200:totYrs,:,:])
meanPburn      = AvgGridsData[56,:,:] # np.mean(AllGridsData[53,200:totYrs,:,:])

#%% --------------------- Plot ---------------------
PFTID = ['Woody','Grass']
xLong = np.arange(N_Lon)
Longi = np.arange(-179.25+minLon*0.5, -179.25+(maxLon+1)*0.5, 0.5*Resolution)
Y_Lon = 1 #14

meanCAgrass[meanCAgrass>1.0] = 1.0

"""
plt.figure(1) #
plt.clf()

plt.subplot(221)
plt.plot(Longi, meanRain[Y_Lon,:],linewidth=1)
#plt.plot(Longi, meanCAgrass[0,:],linewidth=1)
#plt.legend((PFTID),loc=0,ncol=1)
plt.xlabel('Longitude', fontdict=font)
plt.ylabel('Rainfall (mm m$^{-2}$ yr$^{-1}$)', fontdict=font)
#plt.ylim(0,0.25)
plt.title('(a)',x=0.85, y=0.85, fontdict=font)

plt.subplot(222)
plt.plot(Longi,meanGPP[Y_Lon,:],linewidth=1)
#plt.plot(Longi, meanCAgrass[0,:],linewidth=1)
#plt.legend((PFTID),loc=0,ncol=1)
plt.xlabel('Longitude', fontdict=font)
plt.ylabel('GPP (KgC m$^{-2}$ yr$^{-1}$)', fontdict=font)
#plt.ylim(0,0.25)
plt.title('(b)',x=0.85, y=0.85, fontdict=font)

plt.subplot(223)
plt.plot(Longi, meanPlantC[Y_Lon,:],linewidth=1)
#plt.plot(Longi, meanCAgrass[0,:],linewidth=1)
#plt.legend((PFTID),loc=0,ncol=1)
plt.xlabel('Longitude', fontdict=font)
plt.ylabel('PlantC (KgC m$^{-2}$)', fontdict=font)
#plt.ylim(0,0.25)
plt.title('(c)',x=0.85, y=0.85, fontdict=font)

plt.subplot(224)
plt.plot(Longi, meanCAtree[Y_Lon,:],Longi,meanCAgrass[Y_Lon,:],linewidth=1)
#plt.plot(Longi, meanCAgrass[0,:],linewidth=1)
plt.legend((PFTID),loc=0,ncol=1)
plt.xlabel('Longitude', fontdict=font)
#plt.ylabel('GPP\n (KgC m$^{-2}$ yr$^{-1}$)', fontdict=font)
plt.ylabel('Acrown (m$^{2}$ m$^{-2}$)', fontdict=font)
#plt.ylim(0,0.25)
plt.title('(d)',x=0.85, y=0.85, fontdict=font)
#plt.tick_params(axis='x', labelbottom=False)
"""

fig = plt.figure(2)
plt.clf()
plt.subplot(231)
plt.title('Yearly GPP \n (KgC m$^{-2}$ yr$^{-1}$)')
#plt.imshow(fracBA_EP[:,:,1],interpolation='none',extent=[8,24,0.14,0.04],aspect=150)
plt.imshow(np.flipud(meanGPP),interpolation='none') #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)

#cax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#cax.get_xaxis().set_visible(False)
#cax.get_yaxis().set_visible(False)
#cax.set_frame_on(False)
plt.colorbar(orientation='horizontal')

plt.subplot(232)
plt.title('Plant Carbon \n(KgC m$^{-2}$)')
#plt.imshow(fracBA_EP[:,:,1],interpolation='none',extent=[8,24,0.14,0.04],aspect=150)
plt.imshow(np.flipud(meanPlantC),interpolation='none') #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)
plt.colorbar(orientation='horizontal')

plt.subplot(233)
plt.title('Fraction of woody CA')
#plt.imshow(fracBA_EP[:,:,1],interpolation='none',extent=[8,24,0.14,0.04],aspect=150)
plt.imshow(np.flipud(meanCAtree),interpolation='none') #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)

#cax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#cax.get_xaxis().set_visible(False)
#cax.get_yaxis().set_visible(False)
#cax.set_frame_on(False)
plt.colorbar(orientation='horizontal')

plt.subplot(234)
plt.title('Fraction of grass CA')
#plt.imshow(fracBA_EP[:,:,1],interpolation='none',extent=[8,24,0.14,0.04],aspect=150)
plt.imshow(np.flipud(meanCAgrass),interpolation='none') #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)
plt.colorbar(orientation='horizontal')

plt.subplot(235)
plt.title('Fire Risk')
plt.imshow(np.flipud(meanFrisk),interpolation='none') #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)

#cax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#cax.get_xaxis().set_visible(False)
#cax.get_yaxis().set_visible(False)
#cax.set_frame_on(False)
plt.colorbar(orientation='horizontal')

plt.subplot(236)
plt.title('Fire probability')
plt.imshow(np.flipud(meanPburn),interpolation='none') #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)
plt.colorbar(orientation='horizontal')

# PFT CA coverage
fig = plt.figure(3)
plt.clf()
plt.subplot(241)
plt.title('C4 grass')
#plt.imshow(fracBA_EP[:,:,1],interpolation='none',extent=[8,24,0.14,0.04],aspect=150)
plt.imshow(np.flipud(meanPFTCA[0,:,:]),interpolation='none',vmin=0, vmax=1.0) #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)
plt.colorbar(orientation='horizontal')

plt.subplot(242)
plt.title('C3 grass')
#plt.imshow(fracBA_EP[:,:,1],interpolation='none',extent=[8,24,0.14,0.04],aspect=150)
plt.imshow(np.flipud(meanPFTCA[1,:,:]),interpolation='none',vmin=0, vmax=1.0) #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)
plt.colorbar(orientation='horizontal')

plt.subplot(243)
plt.title('Tropical evergreen')
#plt.imshow(fracBA_EP[:,:,1],interpolation='none',extent=[8,24,0.14,0.04],aspect=150)
plt.imshow(np.flipud(meanPFTCA[2,:,:]),interpolation='none') #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)
plt.colorbar(orientation='horizontal')

plt.subplot(244)
plt.title('Topical drought-deciduous')
#plt.imshow(fracBA_EP[:,:,1],interpolation='none',extent=[8,24,0.14,0.04],aspect=150)
plt.imshow(np.flipud(meanPFTCA[3,:,:]),interpolation='none') #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)
plt.colorbar(orientation='horizontal')

plt.subplot(245)
plt.title('Temperate/Boreal evergreen')
#plt.imshow(fracBA_EP[:,:,1],interpolation='none',extent=[8,24,0.14,0.04],aspect=150)
plt.imshow(np.flipud(meanPFTCA[4,:,:]),interpolation='none') #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)
plt.colorbar(orientation='horizontal')

plt.subplot(246)
plt.title('Cold-deciduous')
#plt.imshow(fracBA_EP[:,:,1],interpolation='none',extent=[8,24,0.14,0.04],aspect=150)
plt.imshow(np.flipud(meanPFTCA[5,:,:]),interpolation='none',cmap='viridis', vmin=0, vmax=1.0) #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)
plt.colorbar(orientation='horizontal')

plt.subplot(247)
plt.title('Desert Shrub')
#plt.imshow(fracBA_EP[:,:,1],interpolation='none',extent=[8,24,0.14,0.04],aspect=150)
plt.imshow(np.flipud(meanPFTCA[7,:,:]),interpolation='none',cmap='viridis', vmin=0, vmax=0.6) #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)
plt.colorbar(orientation='horizontal')

plt.subplot(248)
plt.title('N-fixer')
#plt.imshow(fracBA_EP[:,:,1],interpolation='none',extent=[8,24,0.14,0.04],aspect=150)
plt.imshow(np.flipud(meanPFTCA[6,:,:]),interpolation='none',cmap='viridis', vmin=0, vmax=1.0) #,extent=[4,40,0.14,0.04],aspect=300)
#ax.set_aspect('equal')
plt.ylabel('Latitude', fontdict=font)
plt.xlabel('Longitude', fontdict=font)
plt.colorbar(orientation='horizontal')
