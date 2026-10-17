#!/usr/bin/env python3
"""Convert the gzipped classic-netCDF CRU-JRA forcing files used by CRU_Interp.F90
into internally compressed netCDF-4 files (deflate + shuffle), chunked as
(all times, 30 lat, 30 lon) so that a lon/lat sub-block is read without
decompressing the whole file. The Fortran code reads the result directly with
nf90_open; since <name>.nc exists, unzip_gzip_file() skips the gunzip step.

The output mirrors the input layout and the input files are never modified:
  <src>/crujra3/<ncversion><var>.<yr>.365d.noc.nc.gz  -> <dst>/crujra3/...nc
  <src>/Radiation/Processed/tswrf_v15_<yr>.nc.gz      -> <dst>/Radiation/Processed/...nc
SRC_DIR and DST_DIR default to DEFAULT_SRC and DEFAULT_DST (set near the top of this file).
Existing outputs are skipped, so the script can be stopped and restarted.
Then point ncfilepath at <dst> (or move the .nc files next to the .gz files).

Usage:
  python3 ConvertNCtoNC4.py [SRC_DIR] [DST_DIR] [--ncversion crujra.v4.10.5d.]
          [--yr0 1901] [--yr1 2025] [--jobs 6] [--level 1] [--chunk 30]
          [--vars tmp pre spfh pres ugrd vgrd tswrf] [--tmpdir /path/on/fast/disk]
Each worker holds one variable (about 1.5 GB) in memory, so jobs x 1.5 GB of RAM is
needed, plus the same size of temporary space on --tmpdir for the unzipped file.
"""
import argparse
import gzip
import os
import shutil
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

import netCDF4 as nc

ncversion   = 'crujra.v4.10.5d.'
DEFAULT_SRC = '/media/eweng/Expansion/Data/TrendyForcing2025'  # folder with the .nc.gz files
DEFAULT_DST = '/media/eweng/HD2/weng/Data/TrendyForcingUnzipped'  # target directory for the converted .nc files


def to_nc4(src, dst, level, chunk):
    """Copy a classic netCDF file to netCDF-4 with deflate, shuffle and lon/lat chunking."""
    with nc.Dataset(src) as a, nc.Dataset(dst, 'w', format='NETCDF4') as b:
        b.setncatts({k: a.getncattr(k) for k in a.ncattrs()})
        for name, d in a.dimensions.items():
            b.createDimension(name, None if d.isunlimited() else len(d))
        for name, v in a.variables.items():
            kw = {}
            if '_FillValue' in v.ncattrs():
                kw['fill_value'] = v.getncattr('_FillValue')
            if v.ndim == 3:  # (time, lat, lon): chunk = all times x chunk x chunk
                kw.update(zlib=True, complevel=level, shuffle=True,
                          chunksizes=(v.shape[0], min(chunk, v.shape[1]), min(chunk, v.shape[2])))
            elif v.ndim > 0:
                kw.update(zlib=True, complevel=level, shuffle=True)
            o = b.createVariable(name, v.dtype, v.dimensions, **kw)
            o.setncatts({k: v.getncattr(k) for k in v.ncattrs() if k != '_FillValue'})
            v.set_auto_maskandscale(False)  # copy raw values, no masking or scaling
            o.set_auto_maskandscale(False)
            o[...] = v[...]


def convert(job):
    gz, out, tmpdir, level, chunk = job
    if os.path.exists(out):
        return out, 'skipped (exists)'
    os.makedirs(os.path.dirname(out), exist_ok=True)
    part = out + '.part'
    fd, raw = tempfile.mkstemp(suffix='.nc', dir=tmpdir)
    os.close(fd)
    try:
        with gzip.open(gz, 'rb') as fi, open(raw, 'wb') as fo:
            shutil.copyfileobj(fi, fo, 16 * 1024 * 1024)
        to_nc4(raw, part, level, chunk)
        os.replace(part, out)  # atomic: a finished name is always a complete file
    finally:
        for f in (raw, part):
            if os.path.exists(f):
                os.remove(f)
    return out, 'done %.0f MB' % (os.path.getsize(out) / 1e6)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('src', nargs='?', default=DEFAULT_SRC)
    p.add_argument('dst', nargs='?', default=DEFAULT_DST)
    p.add_argument('--ncversion', default=ncversion)
    p.add_argument('--yr0', type=int, default=1901)
    p.add_argument('--yr1', type=int, default=2025)
    p.add_argument('--vars', nargs='+', default=['tmp', 'pre', 'spfh', 'pres', 'ugrd', 'vgrd', 'tswrf'])
    p.add_argument('--jobs', type=int, default=6)
    p.add_argument('--level', type=int, default=1, help='deflate level 1-9')
    p.add_argument('--chunk', type=int, default=30, help='lat/lon chunk size')
    p.add_argument('--tmpdir', default=None, help='where to put the unzipped file (fast local disk)')
    a = p.parse_args()

    jobs, missing = [], []
    for var in a.vars:
        for yr in range(a.yr0, a.yr1 + 1):
            if var == 'tswrf':
                rel = 'Radiation/Processed/tswrf_v15_%d.nc' % yr
            else:
                rel = 'crujra3/%s%s.%d.365d.noc.nc' % (a.ncversion, var, yr)
            gz = os.path.join(a.src, rel + '.gz')
            if not os.path.exists(gz):
                missing.append(gz)
                continue
            jobs.append((gz, os.path.join(a.dst, rel), a.tmpdir, a.level, a.chunk))
    if missing:
        print('Missing %d input files, e.g. %s' % (len(missing), missing[0]))
        sys.exit(1)

    print('Converting %d files with %d workers' % (len(jobs), a.jobs), flush=True)
    failed = 0
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(convert, j): j for j in jobs}
        for n, f in enumerate(as_completed(futs), 1):
            try:
                out, msg = f.result()
                print('[%d/%d] %s: %s' % (n, len(jobs), os.path.basename(out), msg), flush=True)
            except Exception as e:  # keep going; report at the end
                failed += 1
                print('[%d/%d] FAILED %s: %s' % (n, len(jobs), futs[f][0], e), flush=True)
    if failed:
        print('%d file(s) failed. Rerun to retry them.' % failed)
        sys.exit(1)
    print('All files converted.')


if __name__ == '__main__':
    main()
