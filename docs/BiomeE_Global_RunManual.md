# BiomeE Global — Run Manual

This manual covers how to compile and run BiomeE's Global branch for:
- **Single-site (grid cell) runs** using `runBiomeE.x`
- **Generating pre-interpolated forcing CSV files** using `GlobalDataRun.x`
- **Global or regional runs** using `runMultiBlocks.x`

---

## 1. Repository Layout

```
BiomeE-Global/
├── src/                    # Fortran source files
│   ├── datatypes.F90       # Data types, constants, PFT definitions
│   ├── model_utils.F90     # Utility functions
│   ├── io_mod.F90          # I/O, namelist reading, output writers
│   ├── netcdf_io.F90       # NetCDF reading (global mode only)
│   ├── soil.F90            # Soil physics and biogeochemistry
│   ├── vegetation.F90      # Plant physiology, growth, demography
│   ├── animal.F90          # Animal functional types (enabled with -DDO_ANIMAL)
│   ├── restart_mod.F90     # Checkpoint/restart support
│   ├── BiomeE.F90          # Top-level model subroutine
│   └── main.F90            # Entry point (program BiomeE)
├── para_files/             # Namelist (.nml) parameter files
│   ├── parameters_ORNL_test.nml        # Single-site example (ORNL, C3/C4)
│   ├── parameters_ORNL_animal_test.nml # Single-site with animal functional types
│   ├── parameters_GlobalData.nml       # Template for GlobalDataRun.x (WriteForcing)
│   ├── parameters_GlobalBlock.nml      # Template for runMultiBlocks.x
│   ├── parameters_Global_test.nml      # Global run test/debug namelist
│   ├── parameters_DroughtPaleo.nml     # Paleo-drought scenario
│   ├── parameters_TropicalHydro.nml    # Tropical vegetation hydraulics test
│   └── parameters_WIEMIP.nml          # WIEMIP land-use/climate scenario
├── input/                  # Forcing data (single-site CSV/TSV)
├── output/                 # Default output directory
├── runBiomeE.x             # Compile-and-run script: single site
├── GlobalDataRun.x         # Compile-and-run script: generate interpolated forcing CSVs
├── runMultiBlocks.x        # Compile-and-run script: global/regional
└── RunBiomeE.sh            # Modular helper build script
```

---

## 2. Dependencies

| Dependency | Notes |
|---|---|
| `gfortran` ≥ 9 | GNU Fortran compiler |
| NetCDF-Fortran | Required for global runs only (`-lnetcdff`) |
| Standard POSIX shell | `sh`/`bash` for the run scripts |

NetCDF headers and libraries are expected at `/usr/local/include` and `/usr/local/lib` (adjustable in the script headers).

---

## 3. Preprocessor Flags

The model uses C preprocessor (`-cpp`) flags to switch features on and off at compile time.

| Flag | Effect |
|---|---|
| `-DGlobalRun` | Enable global/regional mode (netCDF forcing by defaut, grid loop) |
| `-DDO_Climate_VEG` | Initialize PFT fractions from a climate-envelope map |
| `-DUse_InterpolatedData` | Read pre-interpolated hourly grid csv files |
| `-DZippedNCfiles` | Read gzip-compressed netCDF files |
| `-DZip_outputs` | Compress output files with gzip |
| `-DScreenOutput` | Print diagnostics to stdout |
| `-DDO_ANIMAL` | Enable animal functional types |
| `-DHydro_test` | Hydrology-only test mode |

---

## 4. Single-Site Run with `runBiomeE.x`

### 4.1 What it does

Compiles the model **without** `-DGlobalRun`, links **no** NetCDF libraries, then runs the executable against a single namelist file and a local CSV/TSV forcing file.

### 4.2 Quick start

```bash
# Edit the paths at the top of the script if needed, then:
sh runBiomeE.x
```

The script compiles to an executable named `ess` (or similar) and immediately runs it.

### 4.3 Specifying a different namelist

The executable accepts a namelist path as its first argument:

```bash
./ess ./para_files/my_site.nml
```

If no argument is given, it defaults to `./para_files/input.nml`.

### 4.4 Forcing data format

The forcing file is plain text, tab-separated, one row per hour.

**Required columns (13 fields):**

| Column | Variable | Units |
|---|---|---|
| 1 | YEAR | calendar year |
| 2 | DOY | day of year (1–365) |
| 3 | HOUR | hour of day (0–23) |
| 4 | PAR | mol m⁻² day⁻¹ (set 0 for sub-daily) |
| 5 | Swdown | shortwave down radiation (W m⁻²) |
| 6 | TEMP | air temperature (°C) |
| 7 | SoilT | soil temperature (°C) |
| 8 | RH | relative humidity (%) |
| 9 | RAIN | precipitation (mm per timestep) |
| 10 | WIND | wind speed (m s⁻¹) |
| 11 | PRESSURE | atmospheric pressure (Pa) |
| 12 | aCO2_AW | ambient CO₂ (ppm) |
| 13 | amb_co2 | alternate CO₂ field (ppm) |

Place the file in `input/` and set `climfile` in the namelist (see §6).

### 4.5 Typical single-site namelist (`para_files/parameters_ORNL_test.nml`)

Key sections to configure:

```fortran
&initial_state_nml
  filepath_in   = './input/'
  filepath_out  = './output/'
  climfile      = 'ORNL_forcing.txt'    ! forcing file name inside filepath_in
  runID         = 'ORNL_test_'          ! prefix for all output files

  model_run_years = 500                 ! years to simulate
  CO2_c           = 400.0              ! fixed CO2 (ppm) if not from forcing
  Sc_prcp         = 1.0                ! precipitation scaling factor
  Sc_dT           = 0.0               ! temperature offset (°C)

  do_fire          = .True.
  outputdaily      = .False.
  do_restart_write = .False.

  ! Initial vegetation state
  init_cohort_N     = 1               ! number of starting cohorts
  init_cohort_sps   = 3               ! PFT index (0-based, see §7.2)
  init_cohort_Indiv = 0.002           ! individual density (ind m⁻²)
  init_cohort_bsw   = 0.2            ! initial sapwood C (kgC m⁻²)

  ! Initial soil C-N
  init_fast_SOC   = 0.5              ! fast SOM (kgC m⁻²)
  init_slow_SOC   = 12.0             ! slow SOM (kgC m⁻²)
  init_mineralN   = 5.0E-3           ! mineral N (kgN m⁻²)
  N_input         = 2.0E-3           ! annual N deposition (kgN m⁻² yr⁻¹)
/
```

---

## 5. Generating Pre-Interpolated Forcing Data with `GlobalDataRun.x`

### 5.1 Purpose and workflow

Running `runMultiBlocks.x` with `-DUse_InterpolatedData` is substantially faster than
reading the raw CRUJRA netCDF files at runtime, because the model skips the per-grid
temporal interpolation step and reads ready-made hourly CSV files instead.

`GlobalDataRun.x` is the preprocessing step that produces those CSV files. It compiles
the model with `-DGlobalRun -DDO_Climate_VEG -DZip_outputs` (but **without**
`-DUse_InterpolatedData`) and sets `WriteForcing = .True.` in the namelist. Instead of
running the ecological model, it reads the raw CRUJRA netCDF data, interpolates each
grid cell's climate to hourly resolution, and writes the result to a per-grid CSV file.
Those files are then used by `runMultiBlocks.x`.
When `WriteForcing = .True.`, the model processes will be skipped.

**Two-step workflow:**

```
Step 1: GlobalDataRun.x
  raw CRUJRA netCDF  →  hourly forcing CSVs  (one file per grid cell)

Step 2: runMultiBlocks.x  (compiled with -DUse_InterpolatedData)
  hourly forcing CSVs  →  BiomeE model output
```

### 5.2 Quick start

```bash
# Edit the paths inside the script, then:
sh GlobalDataRun.x
```

The script:
1. Compiles `ess_global` with `-DGlobalRun -DDO_Climate_VEG -DZip_outputs`.
2. Splits the globe into **7 longitude bands** (`Lon1`/`Lon2` arrays) and runs them
   **sequentially** (one band at a time, not in parallel).
3. For each band, generates a block-specific namelist from
   `para_files/parameters_GlobalData.nml` by substituting the `LonStart`, `LonEnd`,
   `GlobalVegGridList`, and `TargetDir` placeholders via `sed`.
4. Runs `./ess_global <block_namelist>` for each band in turn.

### 5.3 Key variables to edit before running

Open `GlobalDataRun.x` and adjust:

| Variable | Default | Description |
|---|---|---|
| `Lon1` / `Lon2` arrays | 7 bands covering 1–720 | Longitude index ranges for each band |
| `runTag` | `'InterpolatedData'` | Subdirectory name appended to the base output path |
| `DIRECTORY` | `/media/eweng/HD2/weng/GlobalESSPFTs/$runTag` | Full path where CSV files are written |
| `fp1` | `./para_files/parameters_GlobalData.nml` | Namelist template |

### 5.4 Namelist template (`para_files/parameters_GlobalData.nml`)

The critical setting is `WriteForcing = .True.` in `&global_setting_nml`, which switches
the model into data-writing mode.

```fortran
&initial_state_nml
  filepath_out = 'TargetDir/'       ! replaced by sed with $DIRECTORY
  runID        = 'ESSPT_'
/

&global_setting_nml
  WriteForcing = .True.             ! write hourly forcing CSVs instead of running model

  ncfilepath   = '/media/eweng/HD2/weng/Data/unzippedNC/'
  ncversion    = 'crujra.v2.4.5d.'
  veg_path     = '/media/eweng/HD2/weng/Data/Vegetation/'
  veg_file     = 'pft2011_0.5x0.5.nc'
  int_fpath    = '/media/eweng/HD2/weng/Data/interpolated/'
  int_prefix   = 'crujra.v2.4.5d.'

  GridListFile = 'GlobalVegGridList.csv'   ! replaced by sed per band

  LowerLon = LonStart      ! replaced by sed: start longitude index of this band
  UpperLon = LonEnd        ! replaced by sed: end   longitude index of this band
  LowerLat = 61            ! ≈ 60°S
  UpperLat = 349           ! ≈ 84.75°N

  yr_start = 1990
  yr_end   = 2019

  grid_No1 = 1
  grid_No2 = 90000
  StepLatLon = 1
/
```

Adjust `ncfilepath`, `veg_path`, `int_fpath`, and `yr_start`/`yr_end` to match your
local data layout and the desired climate period.

### 5.5 Output

One CSV file is written per vegetated land grid cell, named by grid coordinates and
placed in `DIRECTORY`. These files contain hourly meteorological forcing at the same
13-column format described in §4.4 and are read directly by `runMultiBlocks.x` when
compiled with `-DUse_InterpolatedData`.

### 5.6 Connecting to `runMultiBlocks.x`

After `GlobalDataRun.x` finishes, point `int_fpath` in
`para_files/parameters_GlobalBlock.nml` to the directory that was set as `DIRECTORY`
here, then compile and run `runMultiBlocks.x` with `-DUse_InterpolatedData` enabled.

---

## 6. Global / Regional Run with `runMultiBlocks.x`

### 6.1 What it does

1. Compiles the model with `-DGlobalRun` (and usually `-DDO_Climate_VEG -DUse_InterpolatedData -DZip_outputs`), linking NetCDF libraries.
2. Divides the global land grid (57,134 cells at 0.5°×0.5°) into **N parallel blocks**.
3. For each block, generates a block-specific namelist by substituting `StartGrid`/`EndGrid` placeholders with actual grid indices via `sed`.
4. Launches each block as a background (`nohup … &`) process.

### 6.2 Prerequisites

- NetCDF-Fortran installed; adjust include/library paths in the script if they differ from `/usr/local`.
- CRUJRA climate data (netCDF, 0.5°×0.5°) accessible at the path set in `ncfilepath`.
- TRENDY vegetation map (`pft2011_0.5x0.5.nc`) accessible at `veg_path`.
- N-deposition files at `ndp_path`.
- Sufficient disk space in the output directory for all grid-cell outputs.

### 6.3 Quick start

```bash
sh runMultiBlocks.x
```

Edit the variables near the top of the script before running:

| Variable | Purpose |
|---|---|
| `MAXGRID` | Total number of land grid cells (default 57,134 for global) |
| `MAXJOBS` | Number of parallel blocks (default 25) |
| Output path | Directory written to each block namelist (replace `/media/eweng/…` with your path) |

### 6.4 Running a regional subset

To restrict the run to a region, edit the `&global_setting_nml` section of the template namelist:

```fortran
&global_setting_nml
  LowerLon = 1     ! westmost 0.5° cell index  (1 = 180°W)
  UpperLon = 720   ! eastmost (720 = 180°E)
  LowerLat = 61    ! southernmost (61 ≈ 60°S)
  UpperLat = 320   ! northernmost (320 ≈ 60°N)
/
```

Longitude index = (longitude + 180) / 0.5 + 1  
Latitude index = (latitude + 90) / 0.5 + 1

Example — Amazon basin only:
```fortran
  LowerLon = 241   !  60°W → index 241
  UpperLon = 360   !  0°   → index 360
  LowerLat = 121   ! 30°S  → index 121
  UpperLat = 181   ! 0°    → index 181
```

### 6.5 Global namelist template (`para_files/parameters_GlobalBlock.nml`)

Key differences from the single-site namelist:

```fortran
&global_setting_nml
  ncfilepath  = '/data/CRUJRA/'
  ncversion   = 'crujra.v2.4.5d.'
  veg_path    = '/data/Vegetation/'
  veg_file    = 'pft2011_0.5x0.5.nc'
  int_fpath   = '/data/Interpolated/'
  int_prefix  = 'crujra.v2.4.5d.'
  yr_start    = 1990
  yr_end      = 2019

  ! These two lines are replaced by runMultiBlocks.x at runtime:
  grid_No1    = StartGrid
  grid_No2    = EndGrid
  StepLatLon  = 1             ! grid step (1 = every cell, 2 = every other)
/

&initial_state_nml
  filepath_out = '/scratch/GlobalRun_2024/'
  runID        = 'Global_'
  model_run_years = 120
/
```

---

## 7. Namelist Reference

### 7.1 `&soil_data_nml`

```fortran
soiltype = 3        ! 1=Sand 2=LoamySand 3=SandyLoam 5=SandyClayLoam
                    ! 6=ClayLoam 7=Clay 9=SiltClayLoam 12=SiltLoam
thksl = 0.1, 0.2, 0.4, 0.8, 1.5   ! layer thicknesses (m), 5 layers
```

### 7.2 `&vegn_parameters_nml` — PFT traits

Each parameter is an array of 8 values, one per PFT (index 0–7).

Default PFT mapping:

| Index | Type |
|---|---|
| 0 | C4 grass |
| 1 | C3 grass |
| 2 | Tropical evergreen tree |
| 3 | Tropical deciduous tree |
| 4 | Temperate/boreal evergreen tree |
| 5 | Temperate deciduous tree |
| 6 | N-fixing shrub |
| 7 | Deciduous shrub |

Key traits:

```fortran
pt        = 1, 0, 0, 0, 0, 0, 0, 0    ! photosynthesis type (1=C4, 0=C3)
lifeform  = 0, 0, 1, 1, 1, 1, 1, 1    ! 0=herbaceous, 1=woody
phenotype = 0, 0, 1, 0, 1, 0, 0, 0    ! 0=deciduous, 1=evergreen
LMA       = ...                         ! leaf mass per area (kgC m⁻²)
LAImax    = ...                         ! maximum LAI
mu0_topL  = ...                         ! baseline mortality rate (yr⁻¹)
rho_wood  = ...                         ! wood density (kgC m⁻³)
gdd_par1, gdd_par2, gdd_par3           ! growing degree-day phenology
IgniteP   = ...                         ! fire flammability (0–1)
R0_Nfix   = ...                         ! N fixation (kgN kgRootC⁻¹ yr⁻¹)
```

### 7.3 `&initial_state_nml` — run control

| Parameter | Description | Default |
|---|---|---|
| `model_run_years` | Simulation length in years | — |
| `post_yrs` | Extra post-processing years | 0 |
| `CO2_c` | Fixed atmospheric CO₂ (ppm) | 370 |
| `Sc_prcp` | Precipitation scaling (1.0 = no change) | 1.0 |
| `Sc_dT` | Temperature offset (°C) | 0.0 |
| `do_fire` | Enable fire | `.True.` |
| `outputdaily` | Write daily output files | `.False.` |
| `do_restart_write` | Write restart checkpoint | `.False.` |
| `do_restart_read` | Start from a restart file | `.False.` |
| `filepath_in` | Forcing data directory | `./input/` |
| `filepath_out` | Output directory | `./output/` |
| `climfile` | Forcing filename (single-site) | — |
| `runID` | Output filename prefix | — |

---

## 8. Output Files

All output files are written to `filepath_out` with names prefixed by `runID`.

| File suffix | Contents | Frequency |
|---|---|---|
| `*_vegn.csv` | Cohort-level state (biomass, DBH, LAI, mortality) | Annual |
| `*_soil.csv` | Soil C-N pools and fluxes | Annual |
| `*_flux.csv` | Carbon and water fluxes (GPP, NPP, ET, runoff) | Annual |
| `*_daily.csv` | Daily diagnostics | Daily (if `outputdaily=.True.`) |
| `*_restart.bin` | Binary restart file | On demand |

In global mode each grid cell writes its own set of files; outputs may be gzip-compressed if compiled with `-DZip_outputs`.

---

## 9. Restarting a Run

1. Set `do_restart_write = .True.` in the namelist for the initial run to produce `*_restart.bin`.
2. For the continuation run, set:
   ```fortran
   do_restart_read  = .True.
   do_restart_write = .True.   ! to keep writing checkpoints
   ```
3. Point `filepath_in` to the directory containing the restart file, or set the full path explicitly if required.

---

## 10. Common Issues

| Problem | Likely cause | Fix |
|---|---|---|
| `cannot open file` error | Wrong `filepath_in` or `climfile` | Check absolute vs. relative paths |
| NetCDF link error at compile time | Library not found | Set `-I` and `-L` to your NetCDF install |
| All-zero output | Forcing file columns misaligned | Verify tab-separated columns match §4.4 |
| Blocks finish instantly with no output | `grid_No1 > grid_No2` | Check `sed` substitution in `runMultiBlocks.x` |
| High mortality / vegetation collapse | Soil N too low | Increase `init_mineralN` or `N_input` |
| Long spin-up needed | Starting from bare ground | Run ~500 years before the analysis period |
