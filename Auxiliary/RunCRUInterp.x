#!/bin/bash
set -e # Stop at the first failed command (e.g. a block that runs out of memory)
# Self-contained program: only src/CRU_Interp.F90 is needed (no BiomeE simulation is run).
#gfortran src/CRU_Interp.F90 -o cru_interp -I/opt/local/include -L/opt/local/lib -lnetcdff || exit 1
gfortran src/CRU_Interp.F90 -o cru_interp -I/usr/local/include -L/usr/local/lib -lnetcdff || exit 1

# -----------------------------------------------------------------------------
# -------------------Setup data blocks----------------------------------------
runTag='Interpolated20261016' # interpolated forcing files are written here

#! Total grids are 56395 when Lat 61~320 and Lon 1~720 at 0.5x0.5 grid
#Lon1=(1   181 251 381 451 541 621)
#Lon2=(180 250 380 450 540 620 720)
Lon1=(190)
Lon2=(210)

# CRU/vegetation data paths (no namelist file is used)
#ncfilepath='/media/eweng/HD2/weng/Data/unzippedNC/'
#veg_path='/media/eweng/HD2/weng/Data/Vegetation/'
#ncfilepath='/Users/eweng/Documents/Data/CRU/zipped/'
ncfilepath='/media/eweng/Expansion/Data/TrendyForcing2025/'
veg_path='/media/eweng/HD2/weng/Data/Vegetation/'
DIRECTORY=$ncfilepath$runTag

#ncversion='crujra.v2.4.5d.'
ncversion='crujra.v4.10.5d.'
veg_file='pft2011_0.5x0.5.nc'
LowerLat=61
UpperLat=349
yr_start=1901
yr_end=2025
Interpolate=T # T: hourly interpolated forcing; F: write 6-hourly data only (interpolate in BiomeE)

# ----------------- Setup output directory path ------------
echo $DIRECTORY
# Check if the directory exists. If not, create it.
if [ ! -d "$DIRECTORY" ]; then
    echo "Directory $DIRECTORY does not exist. Creating it now..."
    mkdir -p "$DIRECTORY"
    if [ $? -eq 0 ]; then
        echo "Directory $DIRECTORY created successfully."
    else
        echo "Failed to create directory $DIRECTORY."
        exit 1 # Exit with an error code if creation fails
    fi
else
    echo "Directory $DIRECTORY already exists."
fi

# ------------------- Model Run ---------------------------------
for iB in "${!Lon2[@]}"; do
  if [ "${Lon1[$iB]}" -lt "${Lon2[$iB]}" ]; then
    runID='Lon_'${Lon2[$iB]}
    echo "Block ${Lon1[$iB]}-${Lon2[$iB]}"
    echo "Model run: " $runID
    # Run model
    ./cru_interp ncfilepath=$ncfilepath ncversion=$ncversion \
        veg_path=$veg_path veg_file=$veg_file \
        filepath_out=$DIRECTORY/ GridListFile=VegList$runID.csv \
        LowerLon=${Lon1[$iB]} UpperLon=${Lon2[$iB]} \
        LowerLat=$LowerLat UpperLat=$UpperLat yr_start=$yr_start yr_end=$yr_end \
        Interpolate=$Interpolate

  fi
done

rm cru_interp
rm *.mod
