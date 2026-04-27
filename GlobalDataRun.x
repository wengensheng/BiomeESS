#!/bin/bash
FSRCS="src/datatypes.F90 \
       src/model_utils.F90 \
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
#CPPFLAGS+=' -DZippedNCfiles'
#CPPFLAGS+=' -DUse_InterpolatedData'
#CPPFLAGS+=' -DHydro_test'
#CPPFLAGS+=' -DSingleTreeTest'
#CPPFLAGS+=' -DScreenOutput'

echo $FSRCS
echo $CPPFLAGS

#gfortran $FSRCS -o ess -I/opt/local/include -L/opt/local/lib -lnetcdff
#gfortran $FSRCS -o ess -I/Users/eweng/MACPORTS/gcc49-python3/include -L/Users/eweng/MACPORTS/gcc49-python3/lib -lnetcdff
#gfortran src/datatypes.F90 src/io_mod.F90 src/soil.F90 src/vegetation.F90 src/BiomeE.F90 src/main.F90 -DHydro_test -o ess
#gfortran -fopenmp $FSRCS $CPPFLAGS -o ess_global -I/opt/local/include -L/opt/local/lib -lnetcdff

gfortran $FSRCS $CPPFLAGS -o ess_global -I/usr/local/include -L/usr/local/lib -lnetcdff

# -----------------------------------------------------------------------------
# -------------------Setup data blocks----------------------------------------
#! Total grids are 56395 when Lat 61~320 and Lon 1~720 at 0.5x0.5 grid
Lon1=(1   181 251 381 451 541 621)
Lon2=(180 250 380 450 540 620 720)
# namelist file (Parameter and model setting file)
fp1='./para_files/parameters_GlobalData.nml'
echo $fp1

# ----------------- Setup output directory path ------------
runTag='InterpolatedData' #'N3gWmu0Low' #'BaseN2gThnG' #'GrassThn' # 'N2g16Hyrs' #'Warming2C' # 'eCO2'
DIRECTORY="/media/eweng/HD2/weng/GlobalESSPFTs/"$runTag
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
    fp2=$DIRECTORY'/parameters_'$runID'.nml'

    echo "Block ${Lon1[$iB]}-${Lon2[$iB]}"
    echo $fp2
    echo "Model run: " $runID
    sed -e "s/LonStart/${Lon1[$iB]}/g" \
        -e "s/LonEnd/${Lon2[$iB]}/g" \
        -e "s/GlobalVegGridList/VegList$runID/g" \
        $fp1 > $fp2

    echo "Run Longitude ${Lon1[$iB]}-${Lon2[$iB]}"

    # Run model
    ./ess_global $fp2

  fi
done


#rm ess_global
rm *.mod
