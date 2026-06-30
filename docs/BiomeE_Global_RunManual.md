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

## 10. Model Process Description

BiomeE is an individual-based, height-structured vegetation demographic model coupled
to soil carbon–nitrogen biogeochemistry and soil water dynamics. The basic simulation
unit is the **plant cohort** — a group of individuals of the same PFT sharing size,
biomass pools, and physiological state. Cohorts are organized within a **vegetation
tile** (one grid cell), and the tile tracks canopy structure, soil pools, and water
balance. Processes operate on three timescales: **hourly** (physiology, soil water),
**daily** (phenology, growth, tissue turnover), and **annual** (demographics, fire,
hydraulic ageing).

---

### 10.1 Plant Physiology

#### Photosynthesis and Stomatal Conductance

Photosynthesis is computed each hour using the **Leuning (1995) coupled
photosynthesis–stomatal conductance model** implemented in `gs_Leuning`. The canopy is
divided into up to 5 crown layers (CLmax = 5). Light attenuation through the canopy
follows Beer–Lambert extinction, where the PAR fraction reaching layer *i* depends on
the accumulated projected crown area and leaf extinction coefficient of all layers above.

For **C3 species**, gross photosynthesis (Ag) is the minimum of the light-limited rate
(using quantum efficiency α and absorbed PAR) and the Rubisco-limited rate (Farquhar
1980, Vm·(ci−Γ)/(ci+Kc·(1+O/Ko))). For **C4 species** (pt=1), alternative Vm and
light-saturation expressions apply. The maximum carboxylation rate Vm is scaled from
the reference value Vmax by an Arrhenius temperature response (activation energy
24,920 J mol⁻¹). Similarly, Michaelis–Menten constants Kc and Ko are temperature
dependent.

Net photosynthesis An = Ag − Rd, where leaf dark respiration Rd scales with the leaf
nitrogen content per unit area (LNA) and the same Arrhenius temperature function,
suppressed at temperatures below 5 °C and above 45 °C. Stomatal conductance gs follows
the Leuning form: gs = m·An/(ci−Γ)/(1+Ds/D0) + b, where Ds is the leaf-to-air vapour
pressure deficit and m_cond (the g1 parameter) sets the plant's water-use strategy.
When potential transpiration demand Ed exceeds water supply ws from the soil, gs and An
are scaled down proportionally. GPP and transpiration per cohort are computed by
integrating over the cohort's total leaf area (Aleaf = LAI × crown area).

#### Transpiration

Potential transpiration demand Ed (mol H₂O m⁻² leaf s⁻¹) is computed from stomatal
and aerodynamic conductances. Realised transpiration is `transp = min(ws, Ed)`, where
`ws` is the water supply per unit leaf area. In the standard (non-hydraulics) mode,
ws is derived directly from soil water availability via `SoilWaterSupply`. In
plant-hydraulics mode (`-DHydro_test`), ws is calculated explicitly from the stem–leaf
water potential gradient and trunk hydraulic conductance.

#### Maintenance Respiration

Maintenance respiration is calculated each hour in `vegn_respiration` for three
components:
- **Leaf respiration** (`r_leaf`): returned directly from the photosynthesis routine
  as the canopy dark respiration rate (An_cl × Aleaf × step_seconds).
- **Stem (sapwood) respiration** (`r_stem`): proportional to cambium area
  (π × DBH × height × 1.2), the sapwood-specific respiration coefficient γ_SW, and the
  Arrhenius temperature factor tf = exp(9000 × (1/298.16 − 1/TairK)).
- **Root respiration** (`r_root`): proportional to fine root nitrogen content (rootN)
  and γ_FR.

All respiration rates are additionally scaled by `fnsc`, a sigmoid function of the NSC
pool relative to a target (3 × (bl_max + br_max)), reducing respiration when carbon
stores are depleted. NPP per cohort per step is GPP − total respiration.

#### Growth Respiration

A construction cost of 50 % of all new biomass is applied during daily growth
(`resg = 0.5 × dBtotal`), drawn from the NSC pool.

---

### 10.2 Phenology

Phenology is updated daily in `vegn_phenology`. **Evergreen** species keep
`status = LEAF_ON` year-round. **Deciduous** species follow a two-threshold scheme:

- **Leaf-on** is triggered when accumulated growing degree days (GDD, base temperature
  T0_gdd) exceed a chilling-modified threshold `gdd_ON = gdd_par1 + gdd_par2 ×
  exp(gdd_par3 × ncd)`, the smoothed daily temperature (tc_pheno, 80/20 exponential
  smoother) is above `tc0_on`, and soil moisture (thetaS) exceeds a species minimum
  (betaON).
- **Leaf-off** is triggered when tc_pheno falls below a dynamic cold threshold
  `Tc_OFF = tc0_off − 5 × exp(−0.05 × (ngd − N0_GD))`, or soil moisture drops below
  betaOFF, after a minimum growing season length.

At leaf-off, `Seasonal_fall` sheds leaves and (for deciduous species) fine roots at a
daily rate (5 % of bl_max and 2.5 % of br_max per day). A fixed fraction (l_fract) of
senesced leaf, root, and grass-stem carbon is retranslocated to NSC; the nitrogen
retranslocation fraction is retransN. The non-retranslocated fraction enters the fine
litter (SOC1/SON1) and coarse litter (SOC2/SON2) pools.

For deciduous grasses, each new growing season the grass cohort is reset: total
plant carbon and nitrogen are redistributed as a new seedling at the density that can
be supported by the available C and N.

---

### 10.3 Plant Growth and Allocation

Growth is computed daily in `vegn_growth` and `fetch_CN_for_growth`.

#### Carbon and Nitrogen mobilisation

Each day, available carbon for growth (Cgrowth) and nitrogen supply (Nsupply) are
drawn from the non-structural pools (NSC and NSN). The draw combines a **demand-pull**
component (leaf and root filling rate LFR_rate × deficit from bl_max and br_max) and a
**surplus-push** component (excess NSC/NSN above a target drained over a residence time
tauNSC). Cgrowth is capped at 2 % of NSC per day; Nsupply is similarly capped.

#### Carbon allocation

Carbon is allocated in the following priority order:

1. **Leaves (dBL) and fine roots (dBR)**: carbon is spent to fill deficits toward
   `bl_max` and `br_max`. The split between leaves and roots is proportional to their
   respective maximum biomass targets, bounded by the maximum fraction f_LFR_max of
   Cgrowth that can go to leaves and roots.
2. **Seeds (dSeed)**: only for canopy-layer cohorts older than AgeRepro. A fixed
   fraction v_seed of the remaining carbon goes to seed production.
3. **Sapwood (dBSW)**: the remainder of Cgrowth after leaves, roots, and seeds.

For grasses, seeds are allocated in all canopy layers.

#### Nitrogen adjustment

If the nitrogen demand for planned leaf, root, and seed growth exceeds Nsupply, all
three are scaled down by the ratio r_N_SD = Nsupply/Ndemand, and the freed carbon
(cc%extraC) is redirected to sapwood growth. This ensures that when nitrogen is limiting,
plants grow thicker stems rather than thin foliage.

#### Nitrogen pool bookkeeping

- Leaf nitrogen: `leafN += dBL / CNleaf0`
- Root nitrogen: `rootN += dBR / CNroot0`
- Seed nitrogen: `seedN += dSeed / CNseed0`
- Sapwood nitrogen: updated from NSN with a fixed fraction (f_N_add × NSN) transferred
  to wood each day, and the balance of N supply after tissue allocation. Any excess
  above the sapwood C:N target (CNwood0) is returned to NSN.

#### Allometry and architecture

After each growth step, `BM2Architecture` updates height, DBH, crown area (Acrown),
and root zone distribution from the total woody biomass (bsw + bHW). Crown area sets
the maximum leaf biomass `bl_max = f_CO2 × LAImax × LMA × Acrown × (1 − f_cGap)`.
Maximum root biomass `br_max` is derived from `bl_max` by the leaf-to-root ratio.
Maximum NSN is set proportional to the nitrogen needed for full leaf and root growth.

#### Tissue turnover

Daily turnover in `vegn_tissue_turnover` sheds leaves at a rate that accelerates with
leaf age (up to 20 % d⁻¹), fine roots at the species-specific rate α_FR/365, and
grass stems at the leaf rate. Retranslocated C and N return to NSC and NSN; the rest
enters litter pools.

---

### 10.4 Crown Organisation and Canopy Layering

Cohorts are sorted annually by height and assigned to discrete **crown layers** (1 =
top, increasing downward) in `vegn_RelayerCohorts`. The layer assignment determines
the light available to each cohort (Beer–Lambert extinction through all layers above)
and feeds back into the mortality, growth, and maximum leaf area calculations.

The LAI within each layer is tracked (`LAI_L`) and accumulated crown area index (CAI)
determines gap fraction and radiation penetration. Crown area per cohort equals
`nindivs × Acrown`. Cohorts in lower layers receive less light, suppressing their
photosynthesis and increasing mortality (via the layer-dependent factor f_L in
`mortality_rate`).

After demographics, similar cohorts (same species, same layer, similar biomass and
density within tolerance diff_S0) are **merged** by `vegn_mergecohorts` to limit the
total number of cohorts and keep the simulation tractable.

---

### 10.5 Demographic Processes

#### Natural Mortality

Mortality rate (yr⁻¹) is calculated in `mortality_rate` as:

```
mu = mu_bg + (1 − mu_bg) × mu_hydro
```

Background mortality `mu_bg` (capped at 0.5 yr⁻¹) has three multiplicative components:
- **Size effect** f_D: a U-shaped function of DBH (high for seedlings and large old
  trees, minimum at intermediate size), parameterised by A_DBH, B_DBH, D0mu.
- **Layer effect** f_L: mortality increases with layer depth (understory suppression),
  scaled by A_un.
- **Seedling effect** f_S: an exponential decline with DBH (very high for small
  seedlings), scaled by A_sd and B_sd.

Hydraulic failure mortality `mu_hydro` is computed as a logistic function of the annual
transpiration supply/demand ratio (w_scale = annualTrsp / totDemand). When w_scale is
low (chronic water stress), mu_hydro is high. The sensitivity is parameterised by the
species-specific threshold W_mu0.

**Carbon starvation** (`vegn_annual_starvation`) kills an entire cohort instantly if
its NSC drops below 0.01 % of bl_max.

Dead tree C and N are partitioned into fine litter (NSC, seeds, fine roots, leaf cell
wall fraction) and coarse litter (sapwood, heartwood, structural leaf) pools in
`plant2soil`.

#### Reproduction

Each year in `vegn_reproduction`, cohorts in the top canopy layer (layer == 1) that
are older than AgeRepro and have accumulated seed carbon above the minimum seedling
mass s0_plant produce a new cohort of seedlings. Seed carbon and nitrogen pooled from
all reproducible parent cohorts of the same PFT are converted to seedling density:
`nindivs = seedC / s0_plant`. Seedling biomass is initialised by `setup_seedling`:
10 % of totC to fine roots, f_iniBSW × totC to sapwood, and the remainder to NSC;
leaves start at zero (LEAF_OFF).

#### Cohort Management

After each annual cycle, zero-density cohorts are removed, remaining cohorts are
re-sorted into layers, similar cohorts are merged, and empty cohorts are deleted. If
all cohorts go extinct, the vegetation is reset to the initial seedling state.

---

### 10.6 Nitrogen Uptake and Fixation

**Nitrogen uptake** (hourly, `vegn_N_uptake`) uses a Michaelis–Menten equation:
total N uptake rate ρ_N_up = ρ_N_up0 × N_roots / (N_roots0 + N_roots), scaled by the
Arrhenius temperature response of soil at tsoil. Total uptake is proportional to
mineralN, and distributed among cohorts in proportion to their root biomass (only for
cohorts with NSN < NSNmax). Mineral N is decremented by the amount absorbed.

**Nitrogen deposition** is added to mineralN each hour proportional to the annual
deposition rate N_input.

**Biological nitrogen fixation** (`vegn_N_fixation`) is active for PFTs with R0_Nfix > 0
(e.g., N-fixing shrubs). Fixation has an obligate component (a minimum fraction of the
potential rate) and a facultative component that uses surplus carbon (extraC). The
carbon cost of fixation is C0_Nfix kgC per kgN fixed, drawn from NSC and charged to
respiration.

---

### 10.7 Soil Biogeochemical Processes

Soil biogeochemistry is computed hourly in `Soil_BGC` (soil.F90) using a five-pool
coupled carbon–nitrogen model:

| Pool | Symbol | C:N | Description |
|---|---|---|---|
| 1 | SOC1 / SON1 | 50 | Fine (metabolic) litter |
| 2 | SOC2 / SON2 | 150 | Coarse (structural) litter |
| 3 | SOC3 / SON3 | 10 | Microbial biomass |
| 4 | SOC4 / SON4 | 15 | Fast SOM |
| 5 | SOC5 / SON5 | 40 | Slow SOM |

**Decomposition** of litter pools 1 and 2 is a first-order process with rates K0SOM(1)
and K0SOM(2), transferring C and N to the fast and slow SOM pools respectively. SOM
pools 3–5 decay at rates K0SOM(3–5) multiplied by the environmental scalar
`A(tsoil, thetaS)` (a joint temperature–moisture response).

**Microbial growth** from decomposition of pools 4 and 5 is limited by the minimum of
the carbon yield (CUE × d_C) and the nitrogen available at microbial C:N = 10. A
fraction (1 − f_M2SOM) of new microbial C returns to the microbial pool (SOC3);
the rest cycles back to the fast and slow pools.

**Net N mineralisation** is the nitrogen released from decomposing pools minus the
nitrogen incorporated into new microbial biomass. Mineralised N is added to the
mineral nitrogen pool (mineralN) and becomes available for plant uptake.

**Nitrogen losses** from the system:
- **Denitrification** (d_Ngas): proportional to mineralN, scaled by the environmental
  scalar A and the denitrification rate K_DeNitr.
- **Mineral N leaching** (d_Nmin): proportional to mineralN and a runoff-scaled loss
  rate K_rf (which is a saturating function of runoff: K_rf = fdsvN × etaN × runoff /
  (fdsvN + etaN × runoff)).
- **Dissolved organic N (DON) leaching** (dN_SOM4, dN_SOM5): a fraction of the
  decomposed N from fast and slow pools, also scaled by K_rf.

An optional **methane module** (`-DDo_CH4`) partitions a fraction of heterotrophic
respiration into CH4 production under anaerobic conditions (high thetaS), with partial
re-oxidation before emission.

---

### 10.8 Soil Water Dynamics

Soil water is tracked in five layers of configurable thickness (thksl) in
`SoilWaterDynamics` (soil.F90). Each hourly step:

1. **Surface evaporation** is calculated using a Penman–Monteith approach with
   aerodynamic (rAero), canopy (rLAI), and soil surface (rSoil) resistances. rSoil
   increases exponentially as the top-layer moisture approaches the wilting point.
   Surface evaporation is deducted from the top soil layer.

2. **Precipitation infiltration** fills each layer from top to bottom up to field
   capacity (FLDCAP). Any water exceeding field capacity in the bottom layer becomes
   runoff.

3. **Drainage** removes a fraction WaterLeakRate of free water per layer per day,
   passing it to the layer below (or to runoff from the bottom layer).

4. **Transpiration** water is removed from soil layers in proportion to the root area
   index per layer (ArootL) and the soil–root conductance. In standard mode
   (`SoilWaterTranspUpdate`), total transpiration from photosynthesis is apportioned
   by root distribution. In hydraulics mode, water uptake per layer is solved from
   the soil–root–stem water potential gradient.

Soil hydraulic properties (matric potential ψ and conductivity K) for each layer are
updated each hour via `SoilWater_psi_K` using the van Genuchten (or similar) functions
parameterised by soil texture.

Soil wetness for the top three layers (thetaS) feeds back to photosynthesis (water
supply for transpiration), phenology (betaON/OFF thresholds), soil decomposition (A
scalar), and fire risk.

---

### 10.9 Plant Hydraulics

When compiled with `-DHydro_test`, the model replaces the simplified water supply
scheme with explicit plant hydraulic states, updated hourly in
`Plant_water_dynamics_linear`:

- Each cohort tracks separate **leaf water content** (W_lf) and **stem (sapwood) water
  content** (W_sw), converted to water potential via exponential pressure–volume curves:
  ψ = ln(W/Wmax) / CR (where CR is the tissue hydraulic capacitance coefficient).
- **Trunk hydraulic conductance** Ktrunk is the sum of conductances across all sapwood
  rings (up to Ysw_max = 210 years), each degraded by a **percent loss of conductivity
  (PLC)** function of stem water potential: PLC = 1/(1 + (ψ/ψ50)^Kexp).
- Water flows from soil to stem base down the water potential gradient through
  layer-by-layer root–soil conductances k_rs(i) = K_soil(i) × ArootL(i). Stem water
  flows to leaves through Ktrunk.
- Xylem embolism is tracked per ring: accumulated hydraulic usage (accH) and
  embolism-induced damage (plcH) reduce the functional area fraction farea of each ring
  following: `farea = 1 − exp(−r_DF × (1 − (accH + plcH) / WTC0))`. This feeds back
  annually into Ktrunk and the hydraulic failure mortality term.

---

### 10.10 Fire Processes

Fire is evaluated once per year in `vegn_fire`:

1. **Environmental fire risk** (Frisk) is a logistic function of the annual
   precipitation-to-PET ratio (P_ET): `Frisk = 1/(1 + exp(A_MI × (P_ET − MI0Fire)))`.
   Drier years have higher Frisk. Alternatively, Frisk can be fixed at a constant value.

2. **Ignition probability** P_Ign combines the flammabilities of grasses and woody
   plants weighted by their fractional crown cover:
   `P_Ign = 1 − (1 − flmb_G × Frisk) × (1 − flmb_W × Frisk)`,
   where flmb_G = max(IgniteP for grasses) × GrassCA and flmb_W = max(IgniteP for trees) × TreeCA.

3. **Fire occurrence** is stochastic: a random number r_Ign is drawn; fire occurs if
   r_Ign < P_Ign.

4. **Fire effects on vegetation**: each cohort's fire-induced mortality rate is
   `mu_fire = mu0fire × p_fire`. For woody plants in a grass fire, p_fire depends on
   grass biomass (severity s_fireG) and an exponential bark-resistance term (r_BK0 × DBH);
   thicker-barked, larger trees survive better. For woody plants in a canopy fire,
   p_fire depends on tree canopy cover.

5. **Carbon and nitrogen partitioning at fire**:
   - **Volatilised to atmosphere** (Cfire/Nfire): 70 % of leaf C/N, 20 % of
     NSC/NSN and woody C/N from dead plants.
   - **To fine litter** (SOC1/SON1, Cfast/Nfast): 30 % of leaf, all fine roots and
     seeds, 80 % of NSC/NSN.
   - **To coarse litter** (SOC2/SON2, Cslow/Nslow): 80 % of sapwood and heartwood.
   - Surface litter is also partially burned: 70 % of SOC1 and 20 % of SOC2
     volatilised.
   - Fire-released N (Nfire from plants and litter) is added to the mineral N pool.

---

### 10.11 Animal Functional Types (optional, `-DDO_ANIMAL`)

When compiled with `-DDO_ANIMAL`, the model supports animal cohorts within each
vegetation tile. Animals are characterised by diet class (herbivore, carnivore,
omnivore), body mass, intake rates, digestibility, and mortality parameters.

- **Herbivore feeding** removes plant carbon from vegetation cohorts in proportion to
  palatability and plant biomass, following a Michaelis–Menten functional response.
- **Carnivore feeding** preys on other animal cohorts.
- **Excretion and carcasses** return C and N to the fast SOM pool (SOC4/SON4).
- **Starvation mortality** increases when intake falls below the maintenance requirement.
- **Reproduction** is annual, proportional to body condition and r_max.

---

### 10.12 Temporal Integration Summary

| Timescale | Processes |
|---|---|
| Hourly | Photosynthesis, stomatal conductance, transpiration, plant respiration, N uptake, N deposition, N fixation, soil BGC decomposition, soil water dynamics, plant hydraulics (if enabled) |
| Daily | Phenology (leaf-on/off), plant growth and allocation, tissue turnover, leaf senescence, grass thinning, age update, soil water potential |
| Annual | Fire, harvest (if enabled), mortality (background + hydraulic failure + starvation), reproduction, cohort relayering, cohort merging, plant hydraulic state update (ring ageing, xylem embolism), animal reproduction and diagnostics |

---

## 11. Common Issues

| Problem | Likely cause | Fix |
|---|---|---|
| `cannot open file` error | Wrong `filepath_in` or `climfile` | Check absolute vs. relative paths |
| NetCDF link error at compile time | Library not found | Set `-I` and `-L` to your NetCDF install |
| All-zero output | Forcing file columns misaligned | Verify tab-separated columns match §4.4 |
| Blocks finish instantly with no output | `grid_No1 > grid_No2` | Check `sed` substitution in `runMultiBlocks.x` |
| High mortality / vegetation collapse | Soil N too low | Increase `init_mineralN` or `N_input` |
| Long spin-up needed | Starting from bare ground | Run ~500 years before the analysis period |
