#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 07 21:52, 2025

Workstation IP: 192.168.0.31

@author: eweng
"""
import gc
import numpy as np
#import netCDF4 as nc4
import  os
import sys
#from datetime import datetime
import gzip
import csv
import matplotlib.pyplot as plt
from cycler import cycler
import cartopy.crs as ccrs
import cartopy.feature as cfeature

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

# Ecosystem file variables
EcoVars = ['CAI',    'LAI',  'GPP', 'Rauto', 'Rh', 'Burned', \
           'Tavg', 'Rain', 'SoilWater','Transp','Evap','Runoff',\
           'plantC', 'soilC', 'plantN','soilN', \
           'NSC',    'SeedC','leafC', 'rootC', 'swC', 'hwC',\
           'NSN',    'SeedN','leafN', 'rootN', 'swN', 'hwN',\
           'fineL',  'strucL', 'McrbC', 'fastSOC', 'slowSOC',\
           'fineN',  'strucN', 'McrbN', 'fastSON', 'slowSON', 'mineralN',\
           'WC1_5',  'WC2_25', 'WC3_50', 'WC4_100', 'WC5_120', \
           'N_fxed', 'N_uptk', 'Nm_SL', 'Nm_FR','N_loss', \
           'treecover', 'grasscover', 'BMgrass','PET', 'Frisk', 'Pburn',\
           'mu','muC','Indv']

EcoUnit = ['m2/m2','m2/m2','KgC m-2 yr-1','KgC m-2 yr-1','KgC m-2 yr-1','KgC m-2 yr-1',\
           'degree C','mm/year','mm','mm/year','mm/year','mm/year',\
           'KgC m-2','KgC m-2','gN m-2','gN m-2',\
           'KgC m-2','KgC m-2','KgC m-2','KgC m-2','KgC m-2','KgC m-2',\
           'gN m-2','gN m-2','gN m-2','gN m-2','gN m-2','gN m-2',\
           'KgC m-2','KgC m-2','KgC m-2','KgC m-2','KgC m-2',\
           'gN m-2','gN m-2','gN m-2','gN m-2','gN m-2','gN m-2',\
           'mm','mm','mm','mm','mm',\
           'gN m-2 yr-1','gN m-2 yr-1','gN m-2 yr-1','gN m-2 yr-1','gN m-2 yr-1',\
           'm2/m2','m2/m2','KgC m-2','mm/year','times/yr','times/yr', \
           'fraction yr-1','KgC m-2 yr-1','individuals/m2']

subfolder = ['Ecosystem','Cohort']

#%% Check the files
# 'eCO2' # 'N2g1123' #  'N3g1121' # 'N3gLowNfx' # 'N3gTr10' # 'Ndps3g'
# 'N4g1128' #  'Warming2C' # '0.5LonLat_N2g1125'
expID = 'BaseN2gThnG' # 'N3gWmu0Low' # 'TmIgnN3g' # 'MI0Fr2N3g' 
#Resolution = 1 #Grids selected

path0 = '/media/eweng/HD2/weng/GlobalESSPFTs/Simulations/GlobalRun_'+expID+'/'
fpath = path0 +'Yearly/'
fpout = path0 +''

N_pfts   = 8 # total PFTs at one site, 4
N_Layers = 3
Npre     = 6
N_gridV  = 55 # 54, added YealryTmp, 11/23/2025
totYrs   = 0 # Will be updated by reading an ecosystem data file

#%% Read in selected grids file
fselectedgrids = '/media/eweng/HD2/weng/GlobalESSPFTs/SelectedGrids.csv'
#if not os.path.exists(fselectedgrids):
#    fselectedgrids = os.path.join(os.path.dirname(__file__), 'SelectedGrids.csv')

selected_grids = np.genfromtxt(fselectedgrids, delimiter=',', names=True, dtype=None, encoding='utf-8')
if selected_grids.size == 0:
    sys.exit(f'No grids found in {fselectedgrids}')

if selected_grids.shape == ():
    LonLatID_values = [int(selected_grids['LonLatID'])]
else:
    LonLatID_values = [int(x) for x in selected_grids['LonLatID']]

ecofiles = [f'ESSPT_{grid_id:06d}_Ecosystem_yearly.csv.gz' for grid_id in LonLatID_values]
cohfiles = [f'ESSPT_{grid_id:06d}_Cohort_yearly.csv.gz' for grid_id in LonLatID_values]

N_files = len(ecofiles)

#%% Find out Resolution, and file data lines and rows
GridID    = np.zeros((N_files), dtype='int')
LonFiles  = np.zeros((N_files), dtype='int')
LatFiles  = np.zeros((N_files), dtype='int')
VldFiles  = np.zeros((N_files), dtype='int')
# Find out iLon and iLat

N_valid = 0
for ifile in range(N_files):
    VldFiles[ifile] = 1
    N_valid = N_valid + 1
    grid_str = ecofiles[ifile][Npre:Npre+6]
    iLonLat = int(grid_str)
    GridID[ifile]  = iLonLat
    #print (iLonLat)
    iLon = int(iLonLat/1000)
    iLat  = iLonLat - iLon*1000
    LonFiles[ifile]  = iLon
    LatFiles[ifile]  = iLat
    if N_valid == 1:
        fname_gz = fpath + ecofiles[ifile]
        try:
            with gzip.open(fname_gz, 'rt', newline='') as f:
                reader = csv.reader(f)
                header = next(reader)  # Read the header row
                #print(f"Header: {header}")
                LandYrV=(np.array(list(csv.reader(f,delimiter=','))))
        except IOError as e:
            print(f"An I/O error occurred: {e}")
            continue # skip this cycle

        f.close()

        rows = len(LandYrV)
        col  = len(LandYrV[1]) - 1 - 2

N_gridV= col
totYrs = rows

N_Lat = int(N0_Lat)
N_Lon = int(N0_Lon)
YR0   = int(totYrs * 0.2) # for calculating averages


#%% Assign data arrays

# For ecosystem data
AvgGridsData = np.zeros((N_gridV + 3, N_Lat,N_Lon)) # include mu and muC

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

# For cohort data
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

# Temporal arrays for selected grids
grid_mu   = np.full((N_files, totYrs, N_pfts), np.nan)
grid_den  = np.full((N_files, totYrs, N_pfts), np.nan)
grid_CA   = np.full((N_files, totYrs, N_pfts), np.nan)
grid_rain = np.full((N_files, totYrs), np.nan)
grid_PET  = np.full((N_files, totYrs), np.nan)
grid_year = np.arange(1, totYrs + 1, dtype=int)

iRain = EcoVars.index('Rain')
iPET  = EcoVars.index('PET')

#%% Read in ecosystem and cohort data
# Opent output files
feco  = open(fpout + "EcoFileNames.txt", "w")
fcoh  = open(fpout + "CohFileNames.txt", "w")
ftmp  = open(fpout + "CA_Tmp.txt", 'w')
fden  = open(fpout + "Den_Tmp.txt", 'w')
fmu   = open(fpout + "Mu_Tmp.txt", 'w')
iGrid = -1 # Count grids
for ifile in range(N_files):
    iLonLat = GridID[ifile]
    m = LonFiles[ifile] - 1
    n = LatFiles[ifile] - 1

    # Open and read ecosystem file
    fname_gz = fpath + ecofiles[ifile]
    try:
        with gzip.open(fname_gz, 'rt', newline='') as f:
            reader = csv.reader(f)
            header = next(reader)  # Read the header row
            #print(f"Header: {header}")
            LandYrV=(np.array(list(csv.reader(f,delimiter=','))))
    except IOError as e:
        print(f"An I/O error occurred: {e}")
        continue # skip this cycle

    feco.write(ecofiles[ifile] + '\n')

    # Land Data dimensions
    if len(LandYrV) < totYrs/10:
        continue
    rows = len(LandYrV)
    col  = len(LandYrV[1]) - 1
    LandYr = LandYrV[0:rows,0:col].astype(np.float64)
    print (ecofiles[ifile],LonFiles[ifile],LatFiles[ifile])
    AvgGridsData[0:N_gridV, n, m] = np.mean(LandYr[YR0:rows,2:2+N_gridV], axis=0)
    grid_rain[ifile, :rows] = LandYr[:, 2 + iRain]
    grid_PET [ifile, :rows] = LandYr[:, 2 + iPET]

    # Open and read Cohort file
    fname_gz = fpath + cohfiles[ifile]
    try:
        with gzip.open(fname_gz, 'rt', newline='') as f:
            reader = csv.reader(f)
            header = next(reader)  # Read the header row
            #print(f"Header: {header}")
            CCYrV=(np.array(list(csv.reader(f,delimiter=','))))

    except IOError as e:
        print(f"An I/O error occurred: {e}")
        continue # skip this cycle

    fcoh.write(cohfiles[ifile] + '\n')

    # Cohort Data dimensions
    rows = len(CCYrV)
    col  = len(CCYrV[1]) - 1
    CCYr = CCYrV[0:rows,0:col].astype(np.float64)
    totCCL = rows
    print (cohfiles[ifile],totCCL)

    # Calculate PFT-level  GPP, NPP, BA, CA, BM, LAI, height
    iGrid = iGrid + 1

    GPP[:,:] = 0.0
    NPP[:,:] = 0.0 #    = np.zeros((totYrs, N_pfts))
    BA[:,:]  = 0.0 #    = np.zeros((totYrs, N_pfts))
    CA[:,:]  = 0.0 #    = np.zeros((totYrs, N_pfts))
    LA[:,:]  = 0.0 #    = np.zeros((totYrs, N_pfts))
    BM[:,:]  = 0.0 #    = np.zeros((totYrs, N_pfts))
    HT[:,:]  = 0.0 #    = np.zeros((totYrs, N_pfts))
    den[:,:] = 0.0
    mu[:,:]  = 0.0 #    Mortality rate
    muC[:,:] = 0.0 #    Mortality C flux
    for i in range(totCCL-1):
        iYr  = int(CCYr[i,1])-1
        iPFT = int(CCYr[i,4]) # -1
        iLayer = min(2,int(CCYr[i,5])-1)
        # Compute cohort variables
        GPP[iYr,iPFT] = GPP[iYr,iPFT] + CCYr[i,6]*CCYr[i,22]/10000
        NPP[iYr,iPFT] = NPP[iYr,iPFT] + CCYr[i,6]*CCYr[i,23]/10000
        BA[iYr,iPFT]  = BA[iYr,iPFT]  + CCYr[i,6]*PI*0.25*CCYr[i,11]**2
        LA[iYr,iPFT]  = LA[iYr,iPFT]  + CCYr[i,6]*CCYr[i,14]/10000
        BM[iYr,iPFT]  = BM[iYr,iPFT]  + CCYr[i,6]*np.sum(CCYr[i,15:21])/10000
        HT[iYr,iPFT]  = max(HT[iYr,iPFT],CCYr[i,12])
        if iLayer == 0:
            muC[iYr,iPFT] = muC[iYr,iPFT] + CCYr[i,6]*CCYr[i,29]*np.sum(CCYr[i,15:21])/10000
            CA[iYr,iPFT]  = CA[iYr,iPFT]  + CCYr[i,6]*CCYr[i,13]/10000
            den[iYr,iPFT] = den[iYr,iPFT] + CCYr[i,6]
            mu[iYr,iPFT]  = mu[iYr,iPFT]  + CCYr[i,6]*CCYr[i,29] # Density weighted

        # Density weighted
        if (iYr < int(CCYr[i+1,1])-1 or i==totCCL-2):
            #print(iYr,CCYr[i+1,1]-1,i,totCCL-2)
            for j in range(8):
                if den[iYr,j] > 1.0:
                    mu[iYr,j] = mu[iYr,j]/den[iYr,j]
                else:
                    mu[iYr,j] = 0.0

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
    #meanPFTmu [:,n,m] = np.mean(mu [YR0:totYrs,:],axis=0)
    #masked_den = np.where(den[YR0:totYrs, :] >= 1e-4, den[YR0:totYrs, :], np.nan)
    #meanPFTden[:, n, m] = np.nanmean(masked_den, axis=0)
    #masked_mu  = np.where(mu [YR0:totYrs, :] >= 1e-4, mu [YR0:totYrs, :], np.nan)
    #meanPFTmu [:, n, m] = np.nanmean(masked_mu, axis=0)
    for k in range(8):
        if meanPFTden[k,n,m] > 1.0:
            totMu = 0.0
            totDn = 0.0
            for j in range(YR0,totYrs):
                totMu = totMu + den[j, k] * mu[j, k]
                totDn = totDn + den[j, k]
            meanPFTmu [k, n, m] = totMu/totDn

    grid_mu [ifile, :, :] = mu
    grid_den[ifile, :, :] = den
    grid_CA [ifile, :, :] = CA

    # --------------- Write out Grid temporal files ---------------
    # Write temporal CA of the PFTs
    for i in range(N_pfts):
        formatted_row = [f"{num:.2f}" for num in CA[:,i]]
        ftmp.write(",".join(formatted_row) + "\n")

    # Write temporal mu of the PFTs
    for i in range(N_pfts):
        formatted_row = [f"{num:.2e}" for num in mu[:,i]]
        fmu.write(",".join(formatted_row) + "\n")

    # Write temporal density of the PFTs
    for i in range(N_pfts):
        formatted_row = [f"{num:.2f}" for num in den[:,i]]
        fden.write(",".join(formatted_row) + "\n")

    # Remove variables and release memory
    del LandYrV, LandYr, CCYr, CCYrV
    gc.collect()

feco.close()
fcoh.close()
ftmp.close()
fmu.close()
fden.close()

#%% Plot temporal mortality dynamics and relationships with rain/PET
top_n_major_pfts = 2
major_ca_threshold = 0.01

plot_dir = os.path.join(fpout, 'SelectedGridPlots')
os.makedirs(plot_dir, exist_ok=True)

for ifile in range(N_files):
    iLonLat = GridID[ifile]
    m = LonFiles[ifile] - 1 
    n = LatFiles[ifile] - 1

    mean_ca = meanPFTCA[:, n, m]
    order = np.argsort(mean_ca)[::-1]
    major_pfts = [ipft for ipft in order if mean_ca[ipft] > major_ca_threshold][:top_n_major_pfts]

    if len(major_pfts) == 0:
        major_pfts = [int(order[0])]

    years = grid_year.copy()
    rain = grid_rain[ifile, :]
    pet  = grid_PET [ifile, :]

    valid_eco = np.isfinite(rain) & np.isfinite(pet)

    fig = plt.figure(figsize=(12, 11))
    plt.clf()

    ax1 = plt.subplot(311)
    for ipft in major_pfts:
        y = grid_mu[ifile, :, ipft]
        valid = np.isfinite(y)
        ax1.plot(years[valid], y[valid], linewidth=1.5,
                 label=f'{PFT_ID[ipft]} ({PFTnames[ipft]})')
    ax1.set_title(f'Grid {iLonLat:06d}: mortality rate (mu) of major PFTs')
    ax1.set_xlabel('Year')
    ax1.set_ylabel('mu (yr$^{-1}$)')
    ax1.legend(loc='best', fontsize=8)
    ax1.grid(True, alpha=0.3)

    ax2 = plt.subplot(312)
    corr_lines = []
    for ipft in major_pfts:
        y = grid_mu[ifile, :, ipft]
        valid = np.isfinite(y) & np.isfinite(rain)
        if np.count_nonzero(valid) > 2:
            ax2.scatter(rain[valid], y[valid], s=16, alpha=0.7,
                        label=PFT_ID[ipft])
            r = np.corrcoef(rain[valid], y[valid])[0, 1]
            if np.isfinite(r):
                corr_lines.append(f'{PFT_ID[ipft]}: r={r:.2f}')
    ax2.set_title('mu vs Rain')
    ax2.set_xlabel('Rain (mm yr$^{-1}$)')
    ax2.set_ylabel('mu (yr$^{-1}$)')
    if len(major_pfts) > 0:
        ax2.legend(loc='best', fontsize=8)
    if corr_lines:
        ax2.text(0.98, 0.98, '\n'.join(corr_lines), transform=ax2.transAxes,
                 ha='right', va='top', fontsize=8,
                 bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    ax2.grid(True, alpha=0.3)

    ax3 = plt.subplot(313)
    corr_lines = []
    for ipft in major_pfts:
        y = grid_mu[ifile, :, ipft]
        valid = np.isfinite(y) & np.isfinite(pet)
        if np.count_nonzero(valid) > 2:
            ax3.scatter(pet[valid], y[valid], s=16, alpha=0.7,
                        label=PFT_ID[ipft])
            r = np.corrcoef(pet[valid], y[valid])[0, 1]
            if np.isfinite(r):
                corr_lines.append(f'{PFT_ID[ipft]}: r={r:.2f}')
    ax3.set_title('mu vs PET')
    ax3.set_xlabel('PET (mm yr$^{-1}$)')
    ax3.set_ylabel('mu (yr$^{-1}$)')
    if len(major_pfts) > 0:
        ax3.legend(loc='best', fontsize=8)
    if corr_lines:
        ax3.text(0.98, 0.98, '\n'.join(corr_lines), transform=ax3.transAxes,
                 ha='right', va='top', fontsize=8,
                 bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(os.path.join(plot_dir, f'Grid_{iLonLat:06d}_mu_rain_pet.png'),
                dpi=200, bbox_inches='tight')
    #plt.close(fig)

print(f'Plots saved to: {plot_dir}')


#%% Plot global map of selected grid positions
map_fig = plt.figure(figsize=(14, 7))
ax = plt.axes(projection=ccrs.PlateCarree())

# Build lon/lat at grid-cell centers
grid_lon = -180.0 + (LonFiles.astype(float) - 0.5) * 0.5
grid_lat =  -90.0 + (LatFiles.astype(float) - 0.5) * 0.5

ax.set_global()
ax.add_feature(cfeature.LAND, facecolor='0.92', edgecolor='none', zorder=0)
ax.add_feature(cfeature.OCEAN, facecolor='white', edgecolor='none', zorder=0)
ax.add_feature(cfeature.COASTLINE, linewidth=0.8, edgecolor='0.35', zorder=1)
ax.add_feature(cfeature.BORDERS, linewidth=0.4, edgecolor='0.55', zorder=1)
ax.gridlines(draw_labels=True, linewidth=0.3, color='0.75', alpha=1.0, linestyle='-')

ax.scatter(grid_lon, grid_lat, s=60, marker='o', color='crimson',
           alpha=0.9, transform=ccrs.PlateCarree(), zorder=3)

for ifile in range(N_files):
    ax.text(grid_lon[ifile] + 1.5, grid_lat[ifile] + 1.5,
            f'{GridID[ifile]:06d}', fontsize=8,
            transform=ccrs.PlateCarree(), zorder=4)

ax.set_title('Locations of selected grids')
plt.tight_layout()
map_file = os.path.join(plot_dir, 'SelectedGrids_global_map.png')
map_fig.savefig(map_file, dpi=200, bbox_inches='tight')
#plt.close(map_fig)

print(f'Global map saved to: {map_file}')

