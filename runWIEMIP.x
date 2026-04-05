#!/bin/bash
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
CPPFLAGS+=' -DZip_outputs'
CPPFLAGS+=' -DZippedNCfiles'
CPPFLAGS+=' -DRead_Ndps_files'
#CPPFLAGS+=' -DUse_InterpolatedData'
#CPPFLAGS+=' -DWIEMIP_PFT_setting'
#CPPFLAGS+=' -DHydro_test'
#CPPFLAGS+=' -DSingleTreeTest'
#CPPFLAGS+=' -DScreenOutput'

echo $FSRCS
echo $CPPFLAGS

gfortran $FSRCS $CPPFLAGS -o ess_global -I/usr/local/include -L/usr/local/lib -lnetcdff
rm *.mod

# -----------------------------------------------------------------------------
# ----------------- Setup output directory path ------------
GridRS='1' # Grid resolution, 1 for grid by grid (1x1), 2 for skipping one for each lon and lat (2x2)
runTag='WIEMIP'
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

START_VAL=1
if [[ "$GridRS" == "1" ]]; then
    # For 0.5x0.5 grid (when GridRS='1')
    END_VAL=55001
    INCREMENT=2200 # 2500 # Optional step value
else
    # For 1x1 grid (when GridRS='2')
    END_VAL=14001
    INCREMENT=500 # 800 # Optional step value
fi
Grid1=($(seq $START_VAL $INCREMENT $END_VAL))
for (( i=0; i<${#Grid1[@]}; i++ )); do
    Grid2[$i]=$(( Grid1[$i] + INCREMENT - 1 ))
done
echo "Grid1: ${Grid1[@]}"
echo "Grid2: ${Grid2[@]}"

# Access individual elements
#echo "First element: ${Grid1[0]}"

# -----------------------------------------------------------------------------
# ----------- Generate runscript and parameter files for each block ----
fr1='./BiomeEBlockRun.x'
fp1='./para_files/parameters_WIEMIP.nml'
echo "Runscript: " $fr1
echo "nml file: " $fp1

for iB in "${!Grid2[@]}"; do
  runID='Block'$iB
  nmlID='Block'${Grid2[$iB]}
  fr2='./run'$runID'.x'
  fp2=$DIRECTORY'/parameters_'$nmlID'.nml'

  echo "Block ${Grid1[$iB]}-${Grid2[$iB]}"

  if [ "${Grid1[$iB]}" -lt "${Grid2[$iB]}" ]; then
    # Setup the runscript file
    echo "Runscript: " $fr2
    sed -e "s/Grid2/${Grid2[$iB]}/g" \
           $fr1 > $fr2

    # namelist file (Parameter and model setting file)
    echo "nml file ID: " $nmlID
    sed -e "s/StartGrid/${Grid1[$iB]}/g" \
        -e "s/EndGrid/${Grid2[$iB]}/g" \
        -e "s/ResGrid/$GridRS/g" \
        -e "s#TargetDir#$DIRECTORY#g" \
           $fp1 > $fp2

    # Run the model with newly generated runscript and parameter file
    echo "Run block ${Grid1[$iB]}-${Grid2[$iB]}"
    # BiomeE is hardwired to read "input.nml" from "./para_files/"
    cat $fp2 > ./para_files/input.nml
    chmod u+x $fr2
    nohup ./$fr2 > $runTag'_'$runID'.out' 2>&1 &

    # wait 10 seconds before removing "input.nml" and starting
    # the next model run so that the model has time to get started
    sleep 10
    rm $fr2
    rm ./para_files/input.nml
  fi

done
