#!/usr/bin/env python3
"""Python version of src/CRU_Interp.F90 + RunCRUInterp.x.

Reads the CRU-JRA 6-hourly netCDF files, the vegetation map, and writes one forcing
file per vegetated grid cell:
  INTERPOLATE = True : <NC_VERSION><GridID>_forcing.csv.gz    (hourly, 7 columns)
  INTERPOLATE = False: <NC_VERSION><GridID>_forcing6h.csv.gz  (6-hourly, 8 columns, no
                       interpolation and no time shifts; see CRU_Interpolation in netcdf_io.F90)
plus one grid list file (GRID_LIST) with the grids of all longitude blocks.

Set the options in the block below and run:  python3 CRU_Interp.py
The input .nc files may be plain (e.g. converted with ConvertNCtoNC4.py) or .nc.gz
(unzipped into TMP_DIR for the read and then removed).
"""
import gzip
import math
import os
import shutil
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

import netCDF4 as nc
import numpy as np

# =============================== SETTINGS ====================================
INTERPOLATE   = False # True      # True: hourly interpolated files; False: 6-hourly data only
NC_PATH       = '/media/eweng/HD2/weng/Data/TrendyForcingUnzipped/'  # has crujra3/ and Radiation/Processed/
VEG_PATH      = '/media/eweng/HD2/weng/Data/Vegetation/'
VEG_FILE      = 'pft2011_0.5x0.5.nc'
NC_VERSION    = 'crujra.v4.10.5d.'
OUT_DIR       = '/media/eweng/HD2/weng/Data/Interpolated20261016/'   # forcing files and the grid list go here
GRID_LIST     = 'GlobalVegGridList.csv'   # one grid list file (in OUT_DIR) for all blocks
TMP_DIR       = None      # where .nc.gz files are unzipped (fast local disk); None: system default
YR_START      = 1901
YR_END        = 2025
LON_BLOCKS    = [(1, 720)]     # (first, last) longitude index of each block, run one after another
# The memory is set by TILE_LONS, not by the block width, so one block can cover the globe.
LOWER_LAT     = 61        # 61 = -59.75 deg
UPPER_LAT     = 349       # 349 = 84.75 deg
STEP_LATLON   = 1         # skip grids: 1 = all; 2 = one per 2x2 grids
N_INPUT       = 0.0       # kgN m-2 yr-1, uniform N input to all grids (closed N: 0)
N_JOBS        = 8         # worker processes for interpolation and file writing
TILE_LONS     = 10        # longitude columns held in memory at once (about 5 GB per 1000 grids)
FLOAT_FMT     = '%12.4E'  # format of each output value (Fortran E12.4)
GZIP_LEVEL    = 6
SKIP_EXISTING = True      # skip grids whose forcing file exists (restart); list lines still written
# ===========================================================================

# Constants (from CRU_Interp.F90)
NLON, NLAT, NTIME = 720, 360, 1460       # CRU grid and steps per year
HOURS_PER_YEAR = 365 * 24
HOURS_NCSTEP = HOURS_PER_YEAR // NTIME   # 6
STEPS_PER_DAY_NC = 24 // HOURS_NCSTEP    # 4
LON0, LAT0, WLON, WLAT = -180.0, -90.0, 360.0 / NLON, 180.0 / NLAT
SHSHIFT = 182 * 24                       # Southern Hemisphere data shift (hours), to July 2nd
SECONDS_PER_DAY = 24.0 * 3600.0
PI = 3.1415926
MOL_AIR, MOL_H2O, SOLARC = 28.96e-3, 18.0e-3, 1361.0
NC_FIELDS = ['tmp', 'pre', 'tswrf', 'spfh', 'pres', 'ugrd', 'vgrd']
VEG_TOTAL = 'TOTAL_VEG'
VEG_IDS = ['SHRUBS_BD', 'SHRUBS_BE', 'SHRUBS_ND', 'SHRUBS_NE', 'TREES_BD', 'TREES_BE',
           'TREES_ND', 'TREES_NE', 'GRASS_MAN', 'GRASS_NAT']

_G = {}  # data shared with the worker processes (inherited through fork)


# ------------------------------- reading -------------------------------------
def nc_file(var, yr):
    if var == 'tswrf':
        return NC_PATH + 'Radiation/Processed/tswrf_v15_%d.nc' % yr
    return NC_PATH + 'crujra3/%s%s.%d.365d.noc.nc' % (NC_VERSION, var, yr)


class OpenNC:
    """Open a netCDF file; if only <file>.gz exists, unzip it to TMP_DIR for the read."""
    def __init__(self, path):
        self.path, self.tmp = path, None

    def __enter__(self):
        if not os.path.exists(self.path):
            if not os.path.exists(self.path + '.gz'):
                raise FileNotFoundError(self.path + '(.gz)')
            fd, self.tmp = tempfile.mkstemp(suffix='.nc', dir=TMP_DIR)
            os.close(fd)
            with gzip.open(self.path + '.gz', 'rb') as fi, open(self.tmp, 'wb') as fo:
                shutil.copyfileobj(fi, fo, 16 * 1024 * 1024)
        self.ds = nc.Dataset(self.tmp or self.path)
        self.ds.set_auto_maskandscale(False)  # raw values, like the Fortran nf90_get_var
        return self.ds

    def __exit__(self, *a):
        self.ds.close()
        if self.tmp:
            os.remove(self.tmp)


def read_tile(lon_list):
    """Select the vegetated grids of the given longitude columns and read their climate."""
    lon0, lon1 = lon_list[0], lon_list[-1]
    ys, xs = slice(LOWER_LAT - 1, UPPER_LAT), slice(lon0 - 1, lon1)
    with nc.Dataset(VEG_PATH + VEG_FILE) as dv:
        dv.set_auto_maskandscale(False)
        vegetated = dv[VEG_TOTAL][ys, xs]
        vegcover = np.stack([dv[v][ys, xs] for v in VEG_IDS], axis=-1)  # (lat, lon, 10)
    with OpenNC(nc_file(NC_FIELDS[0], YR_START)) as ds:
        first = ds[NC_FIELDS[0]][0, ys, xs]  # land where the first step < 9999

    ilon = [i for i in lon_list]
    ilat = list(range(LOWER_LAT, UPPER_LAT + 1, STEP_LATLON))
    cells = [(lo, la) for lo in ilon for la in ilat
             if first[la - LOWER_LAT, lo - lon0] < 9999.0 and vegetated[la - LOWER_LAT, lo - lon0] > 5.0]
    if not cells:
        return None
    iy = np.array([la - LOWER_LAT for lo, la in cells])
    ix = np.array([lo - lon0 for lo, la in cells])
    nyr = YR_END - YR_START + 1
    totl = nyr * NTIME
    clim = np.empty((len(cells), len(NC_FIELDS), totl), dtype=np.float32)
    crutime = np.empty(totl)
    for j, var in enumerate(NC_FIELDS):
        for i in range(nyr):
            with OpenNC(nc_file(var, YR_START + i)) as ds:
                a = ds[var][:, ys, xs][:, iy, ix]   # (NTIME, ngrid)
                if var == 'tswrf':   # mean flux over the 6 hours after the stamp -> J m-2 per 6 hours
                    a = a * np.float32(21600.0)
                clim[:, j, i * NTIME:(i + 1) * NTIME] = a.T
                if j == 0:
                    crutime[i * NTIME:(i + 1) * NTIME] = ds['time'][:]
    return dict(cells=cells, vegcover=vegcover[iy, ix, :], clim=clim, crutime=crutime)


# ------------------------------- solar geometry ------------------------------
def calc_cosz(td, lat_deg):
    """cos(solar zenith angle) as calc_solarzen in CRU_Interp.F90 (td in days, vectorized)."""
    rad = PI / 180.0
    hour = (td - np.floor(td)) * 24.0
    latrad = lat_deg * rad
    delta = np.arcsin(-math.sin(rad * 23.450) * np.cos(2.0 * PI * (td + 10.0) / 365.0))
    cosz = np.sin(latrad) * np.sin(delta) + np.cos(latrad) * np.cos(delta) * np.cos(rad * 15.0 * (hour - 12.0))
    return np.maximum(cosz, 0.0)


_tables = {}


def lat_tables(lat):
    """Per-latitude tables: daily maximum SW (J m-2) by doy, and cosz by (doy, hour)."""
    if lat not in _tables:
        doy = np.arange(1, 366)[:, None]
        td = doy + np.arange(96)[None, :] / 96.0  # 15-minute steps through the day
        swmax = (calc_cosz(td, lat) * SOLARC * SECONDS_PER_DAY / 96.0).sum(axis=1)
        swmax = np.concatenate([[0.0], swmax])     # index by doy (1..365)
        cz = calc_cosz(doy + np.arange(24)[None, :] / 24.0, lat)  # (365, 24)
        _tables[lat] = (swmax, cz)
    return _tables[lat]


def esat(t_c):
    return 610.78 * np.exp(17.27 * t_c / (t_c + 237.3))


# ------------------------------- interpolation -------------------------------
def interpolate_grid(g, lon_idx, lat_idx):
    """Hourly forcing of one grid; the same steps as CRU_Interpolation. Returns (L, 7)."""
    crutime = _G['crutime']
    G = g.astype(np.float64)               # (7, Nlines): tmp, pre, tswrf, spfh, pres, ugrd, vgrd
    nl = G.shape[1]
    longi = LON0 + (lon_idx - 0.5) * WLON
    lati = LAT0 + (lat_idx - 0.5) * WLAT
    totl = nl * HOURS_NCSTEP
    totyr = totl // HOURS_PER_YEAR
    year_end_check = totyr * HOURS_PER_YEAR
    assert year_end_check == totl, 'data length is not a whole number of years'

    # Local hour of the first record, and the day of year of each day's first step
    tm_local = math.fmod(crutime[0] * 24.0, 24.0) + math.trunc(longi / 15.0)
    m_shift = int(tm_local % 24.0)
    doy_day = np.mod(np.trunc(crutime[::STEPS_PER_DAY_NC]).astype(int), 365) + 1

    # Fraction of the maximum daily radiation, same value for the day's 4 steps
    swmax_tab, cz_tab = lat_tables(lati)
    swdaily = G[2].reshape(-1, STEPS_PER_DAY_NC).sum(axis=1)
    fd_day = np.clip(swdaily / (swmax_tab[doy_day] + 0.0001), 0.0, 1.0)
    fd = np.repeat(fd_day, STEPS_PER_DAY_NC)

    frac = np.arange(HOURS_NCSTEP) / float(HOURS_NCSTEP)

    def lin(x):
        """Linear interpolation of 6-hourly x to hourly. The last interval runs to the
        first record and its six values are in reverse order, as in the Fortran code."""
        xn = np.append(x[1:], x[0])
        out = x[:, None] + (xn - x)[:, None] * frac[None, :]
        out[-1, :] = out[-1, ::-1]
        return out.ravel()

    tair = lin(G[0])
    spfh = lin(G[3])
    pres = lin(G[4])
    wind = lin(np.sqrt(G[5] ** 2 + G[6] ** 2))
    fdh = lin(fd)
    rain = np.repeat(G[1] / (6.0 * 3600.0), HOURS_NCSTEP)   # mm/s, constant within 6 hours
    # (the last interval's reversed rain values are all equal, so no change is needed)

    # Shift data to local time (move down by m_shift hours)
    if m_shift >= 1:
        tair, spfh, pres, wind, fdh, rain = (np.roll(a, m_shift) for a in (tair, spfh, pres, wind, fdh, rain))

    # Radiation from the solar constant, zenith angle and SW fraction
    rad = SOLARC * np.tile(cz_tab.ravel(), totyr) * fdh   # W/m2

    # Southern Hemisphere: shift up 182 days (all variables except the time columns)
    if lati < 0.0:
        tair, spfh, pres, wind, fdh, rain, rad = (np.roll(a, -SHSHIFT) for a in (tair, spfh, pres, wind, fdh, rain, rad))

    rh = np.clip(spfh / MOL_H2O * MOL_AIR * pres / esat(tair - 273.16), 0.01, 0.99)
    out = np.empty((totl, 7))
    out[:, 0] = rad
    out[:, 1] = tair
    out[:, 2] = rh
    out[:, 3] = rain
    out[:, 4] = wind
    out[:, 5] = pres
    out[:, 6] = N_INPUT
    return out


# ------------------------------- writing -------------------------------------
def write_csv_gz(path, header, data, nper=8760):
    """Write header + rows of data (floats) as 'v1,v2,...,vN,' lines, gzip-compressed."""
    ncol = data.shape[1]
    row = (FLOAT_FMT + ',') * ncol + '\n'
    tmp = path + '.part'
    with gzip.open(tmp, 'wb', compresslevel=GZIP_LEVEL) as f:
        f.write((' ' + header + '\n').encode())
        for i in range(0, data.shape[0], nper):
            blk = data[i:i + nper]
            f.write(((row * blk.shape[0]) % tuple(blk.ravel().tolist())).encode())
    os.replace(tmp, path)   # a finished name is always a complete file


def process_grid(k):
    lon_idx, lat_idx = _G['cells'][k]
    grid_id = lon_idx * 1000 + lat_idx
    name = '%s%06d_forcing%s.csv.gz' % (NC_VERSION, grid_id, '' if INTERPOLATE else '6h')
    path = OUT_DIR + name
    if SKIP_EXISTING and os.path.exists(path):
        return grid_id
    g = _G['clim'][k]
    if INTERPOLATE:
        write_csv_gz(path, 'Swdown,Tair,RH,RAIN,WIND,PRESSURE,Ndpst',
                     interpolate_grid(g, lon_idx, lat_idx))
    else:  # 6-hourly values as read, plus the N input column
        d = np.empty((g.shape[1], 8))
        d[:, :7] = g.T
        d[:, 7] = N_INPUT
        write_csv_gz(path, 'Tair,RAIN,SWdown6h,SPFH,PRESSURE,WINDU,WINDV,Ndpst', d, nper=4 * 365)
    return grid_id


# ------------------------------- main ----------------------------------------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    t0 = time.time()
    assert 1 <= LOWER_LAT <= UPPER_LAT <= NLAT, 'LOWER_LAT/UPPER_LAT must be within 1..360'
    ctx = get_context('fork')
    n_total = 0
    with open(OUT_DIR + GRID_LIST, 'w') as flist:   # one grid list for all blocks
        for lon1, lon2 in LON_BLOCKS:
            print('Block %d-%d' % (lon1, lon2), flush=True)
            lons = list(range(max(lon1, 1), min(lon2, NLON) + 1, STEP_LATLON))
            n_grids = 0
            for t in range(0, len(lons), TILE_LONS):
                tile = lons[t:t + TILE_LONS]
                tt = time.time()
                d = read_tile(tile)
                if d is None:
                    continue
                print('  lon %d-%d: %d grids read in %.1f min' % (tile[0], tile[-1], len(d['cells']),
                                                                 (time.time() - tt) / 60), flush=True)
                _G.update(cells=d['cells'], clim=d['clim'], crutime=d['crutime'])
                with ProcessPoolExecutor(max_workers=N_JOBS, mp_context=ctx) as ex:
                    ids = ex.map(process_grid, range(len(d['cells'])), chunksize=4)
                    for k, grid_id in enumerate(ids):
                        flist.write('%8d' % grid_id + ''.join(',%.4E' % v for v in d['vegcover'][k])
                                    + ',%.4E\n' % N_INPUT)
                flist.flush()
                n_grids += len(d['cells'])
                print('  lon %d-%d done; %d grids so far; %.1f min elapsed' %
                      (tile[0], tile[-1], n_grids, (time.time() - t0) / 60), flush=True)
                _G.clear()
            n_total += n_grids
            print('Block %d-%d: forcing written for %d grids.' % (lon1, lon2, n_grids), flush=True)
    print('Done: %d grids listed in %s%s' % (n_total, OUT_DIR, GRID_LIST), flush=True)


if __name__ == '__main__':
    main()
