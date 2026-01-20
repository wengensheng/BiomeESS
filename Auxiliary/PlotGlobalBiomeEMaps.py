#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 10/09/2025 

Workstation IP: 192.168.0.31
conda activate spyder-env
conda update spyder

@author: eweng
"""
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
from matplotlib.colors import ListedColormap

#discrete color scheme
MapColor = ListedColormap(['white', 'green', 'blue','red'])

font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 12,
        }

PFTID = ['C4G',  'C3G',  'TrE',  'TrD',  'TmE',  'TmD',  'Nfx',  'DeS']

#%% BiomeE simulated ecosystem and PFT variables
# 'N2g1123' # 'N3g1121' # 'N3g1115' #  'N2g1117' #  'N2g1116' # 
expID = 'BaseN2gThnG' # 'N3gWmu0Low' # 'TmIgnN3g' # 'MI0Fr2N3g' # 'GrassThn' # 'BaseN2g1129'
path1='/Users/eweng/Documents/Projects/GlobalBiomeESS/Simulations/Test-'+expID+'/'
fname1 = 'BiomeE_Simu_' + expID + '.nc'

sf1 = path1 + fname1
#sf2 = path1 + fname3
# Open netcdf files
ncfile1 = nc4.Dataset(sf1) 
#ncfile2 = nc4.Dataset(sf2)

# Load variables from sf1
lat = ncfile1.variables['Latitude'][:]         # latitude, degrees
lon = ncfile1.variables['Longitude'][:]        # longitude, degrees

LAI     = ncfile1.variables['LAI'][:]
GPP     = ncfile1.variables['GPP'][:]
Rauto   = ncfile1.variables['Rauto'][:]
PlantC  = ncfile1.variables['plantC'][:]
SoilC   = ncfile1.variables['soilC'][:]

CAgrass = ncfile1.variables['grasscover'][:]
CAtree  = ncfile1.variables['treecover'][:]
Frisk   = ncfile1.variables['Frisk'][:]
Pburn   = ncfile1.variables['Pburn'][:]
Rain    = ncfile1.variables['Rain'][:]
PET     = ncfile1.variables['PET'][:]

PlantN  = ncfile1.variables['plantN'][:]
SoilN   = ncfile1.variables['soilN'][:]
mineralN = ncfile1.variables['mineralN'][:]
N_fxed   = ncfile1.variables['N_fxed'][:]
N_loss   = ncfile1.variables['N_loss'][:]
N_uptk   = ncfile1.variables['N_uptk'][:]

# Plant function types
C4G = ncfile1.variables['C4G'][:]
C3G = ncfile1.variables['C3G'][:]
TrE = ncfile1.variables['TrE'][:]
TrD = ncfile1.variables['TrD'][:]
TmE = ncfile1.variables['TmE'][:]
TmD = ncfile1.variables['TmD'][:]
Nfx = ncfile1.variables['Nfx'][:]
DeS = ncfile1.variables['DeS'][:]

N_PFT = len(PFTID)
N_lat = len(TrE[:])
N_lon = len(TrE[0][:])
CApft = np.zeros((N_PFT,N_lat,N_lon))
for i in range(N_PFT):
    CApft[i,:,:]=ncfile1.variables[PFTID[i]][:]
    
# Close the netcdf file
ncfile1.close
#ncfile2.close

#%% ============= Data analysis ==================

trees = TrE+TrD+TmE+TmD
grass = C4G + C3G
shrub = DeS
MI = Rain/(PET+0.001)

# Calculate the number of species at each grid
N_PFT = len(PFTID)
N_lat = len(TrE[:])
N_lon = len(TrE[0][:])
pfts = np.zeros((N_PFT,N_lat,N_lon), dtype=int)
nsp = np.zeros((N_lat,N_lon), dtype=int)
thd = 0.2
for i in range(N_PFT):
    pfts[i] = (CApft[i,:,:] > thd).astype(int)

nsp = np.sum(pfts, axis=0)


#%%  Plot
# Discrete intervals for the colorbar
#clevs       = np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

# Color for the ocean mask
ocean_color = np.float64([209,230,241])/255

# Create New Colormap With a Limited Number of entries
nmap=plt.cm.get_cmap(name='viridis',lut=9)



#%% Plot PFT distribution
fig = plt.figure(1,figsize=(20, 15),facecolor="white")
plt.clf()
#fig.suptitle('PFT distribution (crown area, m$^{2}$ m$^{-2}$)') #+run_time)


# 1
ax = fig.add_subplot(2, 3, 1, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
m = ax.pcolormesh(lon, lat, TrE, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=1.0)

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(a)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Tropical evergreen trees')

# 2
ax = fig.add_subplot(2, 3, 2, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
m = ax.pcolormesh(lon, lat, TrD, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=1.0)

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(b)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Tropical drought-deciduous trees')

# 3
ax = fig.add_subplot(2, 3, 3, projection=ccrs.Robinson(central_longitude=0, globe=None))
m = ax.pcolormesh(lon, lat, TmE, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=1.0)
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(c)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Cold-evergreen trees')

# 4
ax = fig.add_subplot(2, 3, 4, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
m = ax.pcolormesh(lon, lat, TmD, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=1.0)

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(d)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Cold-deciduous trees')


# 5
ax = fig.add_subplot(2, 3, 5, projection=ccrs.Robinson(central_longitude=0, globe=None))
m = ax.pcolormesh(lon, lat, grass, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=1.0)
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(e)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Grasses')


# 6
ax = fig.add_subplot(2, 3, 6, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
m = ax.pcolormesh(lon, lat, DeS, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=1.0)

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(f)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Dry shrubs')

plt.tight_layout()
plt.show()

#============== Rainfall and fire risks ========
fig = plt.figure(2,figsize=(20, 15),facecolor="white")
plt.clf()
#fig.suptitle('GPP') #+run_time)

# 1
ax = fig.add_subplot(2, 2, 1, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
m = ax.pcolormesh(lon, lat, Frisk, transform=ccrs.PlateCarree(), \
              cmap='Reds', vmin=0.0,vmax=1.0)

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(a)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Climatic fire risk')

# 2
ax = fig.add_subplot(2, 2, 2, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
m = ax.pcolormesh(lon, lat, Pburn, transform=ccrs.PlateCarree(), \
              cmap='Reds', vmin=0.0,vmax=1.0)

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(b)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Fire probability (in a year)')

# 3
ax = fig.add_subplot(2, 2, 3, projection=ccrs.Robinson(central_longitude=0, globe=None))
m = ax.pcolormesh(lon, lat, CAtree, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=1.0)
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(c)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Woody plant coverage (m$^{2}$ m$^{-2}$)', rotation=0)


# 4
ax = fig.add_subplot(2, 2, 4, projection=ccrs.Robinson(central_longitude=0, globe=None))
m = ax.pcolormesh(lon, lat, CAgrass, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=1.0)
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(d)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Grass coverage (m$^{2}$ m$^{-2}$)', rotation=0)

plt.tight_layout()
plt.show()

#%% Plot cold-deciduous and evergreen
fig = plt.figure(3,figsize=(20, 15),facecolor="white")
plt.clf()
#fig.suptitle('PFT distribution (crown area, m$^{2}$ m$^{-2}$)') #+run_time)

# 1
ax = fig.add_subplot(2, 2, 1, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
m = ax.pcolormesh(lon, lat, TrE, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=1.2)

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(a)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Tropical evergreen trees')

# 2
ax = fig.add_subplot(2, 2, 2, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
m = ax.pcolormesh(lon, lat, TrD, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=1.2)

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(b)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Tropical drought-deciduous trees')

# 3
ax = fig.add_subplot(2, 2, 3, projection=ccrs.Robinson(central_longitude=0, globe=None))
m = ax.pcolormesh(lon, lat, TmE, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=1.2)
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(c)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Cold-evergreen trees')

# 4
ax = fig.add_subplot(2, 2, 4, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
m = ax.pcolormesh(lon, lat, TmD, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=1.2)

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(d)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Cold-deciduous trees')

#============== Arid Semi-arid grasses and shrubs ========
fig = plt.figure(4,figsize=(20, 15),facecolor="white")
plt.clf()
#fig.suptitle('GPP') #+run_time)

# 1
ax = fig.add_subplot(1, 2, 1, projection=ccrs.Robinson(central_longitude=0, globe=None))
m = ax.pcolormesh(lon, lat, CAgrass, transform=ccrs.PlateCarree(), \
              cmap='Greens') #, vmin=0.0,vmax=2.0)
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(a)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Grass coverage (m$^{2}$ m$^{-2}$)', rotation=0)
plt.show()

# 4
ax = fig.add_subplot(1, 2, 2, projection=ccrs.Robinson(central_longitude=0, globe=None))
m = ax.pcolormesh(lon, lat, DeS, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=0.8)
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(b)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Dry shrub coverage (m$^{2}$ m$^{-2}$)', rotation=0)

plt.tight_layout()
plt.show()


#%% Ecosystem properties
fig = plt.figure(5,figsize=(20, 15),facecolor="white")
plt.clf()
#fig.suptitle('GPP') #+run_time)

# 1
ax = fig.add_subplot(2, 2, 1, projection=ccrs.Robinson(central_longitude=0, globe=None))
m = ax.pcolormesh(lon, lat, MI, transform=ccrs.PlateCarree(), \
              cmap='Blues', vmin=0.0,vmax=2.0) # , vmin=0.0,vmax=4000
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(a)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Moisture index (precipitation/PET)', rotation=0)

# 2
ax = fig.add_subplot(2, 2, 2, projection=ccrs.Robinson(central_longitude=0, globe=None))
m = ax.pcolormesh(lon, lat, LAI, transform=ccrs.PlateCarree(), \
              cmap='Greens')
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(b)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('LAI (m$^{2}$ m$^{-2}$)', rotation=0)


# 3
ax = fig.add_subplot(2, 2, 3, projection=ccrs.Robinson(central_longitude=0, globe=None))
m = ax.pcolormesh(lon, lat, GPP, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=3.0)
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(c)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('GPP (kgC m$^{-2}$ yr$^{-1}$)', rotation=0)

# 4
ax = fig.add_subplot(2, 2, 4, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
m = ax.pcolormesh(lon, lat, PlantC, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=30.0)

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(d)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Biomass (kgC m$^{-2}$)')

# 4
#ax = fig.add_subplot(2, 2, 4, projection=ccrs.Robinson(central_longitude=0, globe=None))
##m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
##              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
#m = ax.pcolormesh(lon, lat, SoilC, transform=ccrs.PlateCarree(), \
#              cmap='Greens', vmin=0.0,vmax=40.0)

#ax.coastlines()
#ax.set_global()
#ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
#ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
#ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
#plt.title('(d)', x=0.1, y=0.9)
#cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
#cbar.ax.tick_params(labelsize=12)
#cbar.ax.set_xlabel('Soil carbon (KgC m$^{-2}$)')

plt.tight_layout()
plt.show()

#%% Ecosystem properties - Nitrogen
fig = plt.figure(6,figsize=(20, 15),facecolor="white")
plt.clf()
#fig.suptitle('GPP') #+run_time)

# 1
ax = fig.add_subplot(2, 2, 1, projection=ccrs.Robinson(central_longitude=0, globe=None))
m = ax.pcolormesh(lon, lat, mineralN, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=20.0)
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(a)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Mineral nitrogen (gN m$^{-2}$)', rotation=0)


# 2
ax = fig.add_subplot(2, 2, 2, projection=ccrs.Robinson(central_longitude=0, globe=None))
m = ax.pcolormesh(lon, lat, N_fxed, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=0.5)
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(b)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Fixed nitrogen (gN m$^{-2}$ yr$^{-1}$)', rotation=0)

# 3
ax = fig.add_subplot(2, 2, 3, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
m = ax.pcolormesh(lon, lat, N_uptk, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=20.0)

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(c)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('Plant N update (gN m$^{-2}$ yr$^{-1}$)')

# 4
ax = fig.add_subplot(2, 2, 4, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
m = ax.pcolormesh(lon, lat, N_loss, transform=ccrs.PlateCarree(), \
              cmap='Greens')

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
plt.title('(d)', x=0.1, y=0.9)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_xlabel('N loss rate (gN m$^{-2}$ yr$^{-1}$)')

plt.tight_layout()
plt.show()


#%% ================
fig = plt.figure(21,figsize=(20, 15),facecolor="white")
plt.clf()
#fig.suptitle('Gross Primary Production') #+run_time)

# 1
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0, globe=None))
m = ax.pcolormesh(lon, lat, GPP, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=3.0)
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
#plt.title('(c)', x=0.1, y=0.9)
cbar = plt.colorbar(m, orientation="vertical", fraction=0.04, pad=0.02, shrink=0.8)
cbar.ax.tick_params(labelsize=12)

# Correct way to label a vertical colorbar
cbar.set_label('GPP (kgC m$^{-2}$ yr$^{-1}$)', rotation=90, labelpad=5)



plt.tight_layout()
plt.show()

#% ======== Biomass ============
fig = plt.figure(22,figsize=(20, 15),facecolor="white")
plt.clf()
#fig.suptitle('Biomass') #+run_time)

ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
#m = ax.pcolormesh(lon, lat, Nfx, transform=ccrs.PlateCarree(), \
#              cmap='Greens', vmin=0.0,vmax=1.0)
m = ax.pcolormesh(lon, lat, PlantC, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=30.0)

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
#plt.title('(d)', x=0.1, y=0.9)
# --- COLORBAR ---
cbar = plt.colorbar(m, orientation="vertical", fraction=0.04, pad=0.02, shrink=0.8)
cbar.ax.tick_params(labelsize=12)

# Correct way to label a vertical colorbar
cbar.set_label('Biomass (kgC m$^{-2}$)', rotation=90, labelpad=5)


plt.tight_layout()
plt.show()

#%% ======== PFTs ============
#'C4G',  'C3G',  'TrE',  'TrD',  'TmE',  'TmD',  'Nfx',  'DeS'
# 'trees', 'grass', 'shrub'
fig = plt.figure(23,figsize=(20, 15),facecolor="white")
plt.clf()
#fig.suptitle('Nitrogen fixers') #+run_time)

ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
m = ax.pcolormesh(lon, lat, Nfx, transform=ccrs.PlateCarree(), \
              cmap='Greens', vmin=0.0,vmax=1.0)

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
#plt.title('(d)', x=0.1, y=0.9)
cbar = plt.colorbar(m, orientation="vertical", fraction=0.04, pad=0.02, shrink=0.8)
cbar.ax.tick_params(labelsize=12)
cbar.set_label('CAI (m$^{2}$ m$^{-2}$)', rotation=90, labelpad=5)

plt.tight_layout()
plt.show()

#%% ======== PFT coexistence ============
#'C4G',  'C3G',  'TrE',  'TrD',  'TmE',  'TmD',  'Nfx',  'DeS'
# 'trees', 'grass', 'shrub'
fig = plt.figure(24,figsize=(20, 15),facecolor="white")
plt.clf()
#fig.suptitle('Nitrogen fixers') #+run_time)

ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0, globe=None))
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
m = ax.pcolormesh(lon, lat, nsp, transform=ccrs.PlateCarree(), \
              cmap='Greens')

ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(cartopy.feature.LAND,zorder=1,facecolor='none',edgecolor='k',linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
#plt.title('(d)', x=0.1, y=0.9)
cbar = plt.colorbar(m, orientation="vertical", fraction=0.04, pad=0.02, shrink=0.8)
cbar.ax.tick_params(labelsize=12)
cbar.set_label('N pfts', rotation=90, labelpad=5)

plt.tight_layout()
plt.show()