#!/bin/bash
set -euo pipefail
# -----------------------------------------------------------------------------
# --------- Compile the model --------------------
FSRCS="src/datatypes.F90 \
       src/io_mod.F90 \
       src/netcdf_io.F90 \
       src/soil.F90 \
       src/vegetation.F90 \
       src/BiomeE.F90 \
       src/main.F90"

CPPFLAGS=''
CPPFLAGS+=' -DGlobalRun'
CPPFLAGS+=' -DDO_Climate_VEG'
CPPFLAGS+=' -DUse_InterpolatedData'
#CPPFLAGS+=' -DZippedNCfiles'
#CPPFLAGS+=' -DHydro_test'
#CPPFLAGS+=' -DSingleTreeTest'
#CPPFLAGS+=' -DScreenOutput'
CPPFLAGS+=' -DZip_outputs'

echo $FSRCS
echo $CPPFLAGS

gfortran $FSRCS $CPPFLAGS -o ess_global -I/usr/local/include -L/usr/local/lib -lnetcdff

# Ensure executable built
if [ ! -x "./ess_global" ]; then
  echo "ERROR: ess_global was not created or is not executable." >&2
  exit 1
fi
rm -f *.mod

# -----------------------------------------------------------------------------
# ----------------- Setup output directory path ------------
GridRS='1' # Grid resolution, 1 for grid by grid (1x1), 2 for skipping one for each lon and lat (2x2)
runTag='Test0' #'N3gWmu0Low' #'BaseN2gThnG' #'GrassThn' # 'N2g16Hyrs' #'Warming2C' # 'eCO2'
DIRECTORY="/media/eweng/HD2/weng/GlobalESSPFTs/Simulations/GlobalRun_"$runTag

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

# -----------------------------------------------------------------------------
# ----------- Create block arrays using command substitution with 'seq' ------
#! Total grids are 56395 when Lat 61~320 and Lon 1~720 at 0.5x0.5 grid
#Grid1=(1    3001 6001 9001  12001 15001 18001 21001 24001 27001 30001 33001 36001 39001 42001 45001 48001 51001 54001)
#Grid2=(3000 6000 9000 12000 15000 18000 21000 24000 27000 30000 33000 36000 39000 42000 45000 48000 51000 54000 57500)
#Grid1=(1   801  1601 2401 3201 4001 4801 5601 6401 7201 8001 8801 9601  10401 11201 12001 12801 13601)
#Grid2=(800 1600 2400 3200 4000 4800 5600 6400 7200 8000 8800 9600 10400 11200 12000 12800 13600 14400)

# --- user settings ---
START_VAL=1
MAXGRID=56395
MAXJOBS=25          # number of blocks AND max concurrent jobs (here they match)

# --- derived settings ---
N=$(( MAXGRID - START_VAL + 1 ))
INCREMENT=$(( (N + MAXJOBS - 1) / MAXJOBS ))   # ceil(N/MAXJOBS)

echo "Total grids: $N  MAXJOBS: $MAXJOBS  INCREMENT: $INCREMENT"

# --- build exactly MAXJOBS blocks (or fewer if N < MAXJOBS) ---
Grid1=()
Grid2=()

for ((b=0; b<MAXJOBS; b++)); do
  g1=$(( START_VAL + b*INCREMENT ))
  if (( g1 > MAXGRID )); then
    break
  fi
  g2=$(( g1 + INCREMENT - 1 ))
  if (( g2 > MAXGRID )); then
    g2=$MAXGRID
  fi
  Grid1+=("$g1")
  Grid2+=("$g2")
done
echo "Blocks created: ${#Grid1[@]}"

echo "Grid1: ${Grid1[@]}"
echo "Grid2: ${Grid2[@]}"

# Access individual elements
#echo "First element: ${Grid1[0]}"

# -----------------------------------------------------------------------------
# ----------- Generate runscript and parameter files for each block ----

fp1='./para_files/parameters_GlobalBlock.nml'
echo "nml file: " $fp1

pids=()
for iB in "${!Grid2[@]}"; do
  runID='Block'$iB
  nmlID='Block'${Grid2[$iB]}
  fp2=$DIRECTORY'/parameters_'$nmlID'.nml'

  echo "Block ${Grid1[$iB]}-${Grid2[$iB]}"

  if [ "${Grid1[$iB]}" -le "${Grid2[$iB]}" ]; then

    # namelist file (Parameter and model setting file)
    echo "nml file ID: " $nmlID
    sed -e "s/StartGrid/${Grid1[$iB]}/g" \
        -e "s/EndGrid/${Grid2[$iB]}/g" \
        -e "s/ResGrid/$GridRS/g" \
        -e "s#TargetDir#$DIRECTORY#g" \
           $fp1 > $fp2

    # Run the model with newly generated runscript and parameter file
    echo "Run block ${Grid1[$iB]}-${Grid2[$iB]}"

    # run ess_global
    # nohup ./ess_global $fp2 > $runTag'_'$runID'.out' 2>&1 &
    PROCNAME="${runTag}_${Grid1[$iB]}_${Grid2[$iB]}"
    echo $PROCNAME
    nohup bash -c "exec -a '${PROCNAME}' ./ess_global '${fp2}'" \
      > "${DIRECTORY}/${runTag}_${runID}.out" 2>&1 &
    

    sleep 2
  fi

done

echo "All blocks submitted (MAXJOBS=${MAXJOBS})."