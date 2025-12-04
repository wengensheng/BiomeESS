#!/bin/bash
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

echo $FSRCS
echo $CPPFLAGS

gfortran $FSRCS $CPPFLAGS -o ess_global -I/usr/local/include -L/usr/local/lib -lnetcdff
rm *.mod

# Define the directory path
GridRS='1'
runTag='BaseN2gThnG' #'GrassThn' # 'N2g16Hyrs' #'Warming2C' # 'eCO2'
DIRECTORY="/media/eweng/HD2/weng/GlobalESSPFTs/Simulations/GlobalRun_"$runTag

# Check if the directory does NOT exist
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

#! Total grids are 56395 when Lat 61~320 and Lon 1~720 at 0.5x0.5 grid
#Grid1=(1    3001 6001 9001  12001 15001 18001 21001 24001 27001 30001 33001 36001 39001 42001 45001 48001 51001 54001)
#Grid2=(3000 6000 9000 12000 15000 18000 21000 24000 27000 30000 33000 36000 39000 42000 45000 48000 51000 54000 57500)
#Grid1=(1   801  1601 2401 3201 4001 4801 5601 6401 7201 8001 8801 9601  10401 11201 12001 12801 13601)
#Grid2=(800 1600 2400 3200 4000 4800 5600 6400 7200 8000 8800 9600 10400 11200 12000 12800 13600 14400)
# Create the array using command substitution with 'seq'
START_VAL=1
if [[ "$GridRS" == "1" ]]; then
    # For 0.5x0.5 grid
    END_VAL=55001
    INCREMENT=2500 # Optional step value
else
    # For 1x1 grid
    END_VAL=13601
    INCREMENT=800 # Optional step value 
fi
Grid1=($(seq $START_VAL $INCREMENT $END_VAL))
for (( i=0; i<${#Grid1[@]}; i++ )); do
    Grid2[$i]=$(( Grid1[$i] + INCREMENT - 1 ))
done
echo "Grid1: ${Grid1[@]}"
echo "Grid2: ${Grid2[@]}"

# Access individual elements
#echo "First element: ${Grid1[0]}"

fr1='./BiomeEBlockRun.x'
fp1='./para_files/parameters_GlobalBlock.nml'
echo "Runscript: " $fr1
echo "nml file: " $fp1 

for iB in "${!Grid2[@]}"; do
  runID='Block'$iB
  nmlID='Block'${Grid2[$iB]}
  fr2='./run'$runID'.x'
  #fp2='./para_files/parameters_'$nmlID'.nml'
  fp2=$DIRECTORY'/parameters_'$nmlID'.nml'
  
  echo "Block ${Grid1[$iB]}-${Grid2[$iB]}"
  
  # Setup run script and parameter files
  if [ "${Grid1[$iB]}" -lt "${Grid2[$iB]}" ]; then
    echo "Runscript: " $fr2
    sed -e "s/Grid2/${Grid2[$iB]}/g" $fr1 > $fr2

    # namelist file (Parameter and model setting file)
    echo "nml file ID: " $nmlID
    sed -e "s/StartGrid/${Grid1[$iB]}/g" \
        -e "s/EndGrid/${Grid2[$iB]}/g" \
        -e "s/ResGrid/$GridRS/g" \
        -e "s#TargetDir#$DIRECTORY#g" \
    $fp1 > $fp2
    
    echo "Run block ${Grid1[$iB]}-${Grid2[$iB]}"
    cat $fp2 > ./para_files/input.nml
    chmod u+x $fr2
    nohup ./$fr2 > $runTag'_'$runID'.out' 2>&1 &

    # wait 10 seconds for model being initiated before starting the next model run
    sleep 10
    rm $fr2
    rm ./para_files/input.nml
  fi
  
done
