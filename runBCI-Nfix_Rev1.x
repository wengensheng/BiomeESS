#!/bin/bash

FSRCS="src/datatypes.F90 \
       src/model_utils.F90 \
       src/io_mod.F90 \
       src/restart_mod.F90 \
       src/soil.F90 \
       src/vegetation.F90 \
       src/animal.F90 \
       src/BiomeE.F90 \
       src/main.F90"

CPPFLAGS=''
CPPFLAGS+=' -DZip_outputs'
#CPPFLAGS+=' -DScreenOutput'
#CPPFLAGS+=' -DFixedFireRisk'
#CPPFLAGS+=" -DHydro_test"

#gfortran $FSRCS -o ess -I/opt/local/include -L/opt/local/lib -lnetcdff
#gfortran $FSRCS -o ess -I/Users/eweng/MACPORTS/gcc49-python3/include -L/Users/eweng/MACPORTS/gcc49-python3/lib -lnetcdff
#gfortran src/datatypes.F90 src/io_mod.F90 src/soil.F90 src/vegetation.F90 src/BiomeE.F90 src/main.F90 -DHydro_test -o ess
echo 'gfortran' $FSRCS $CPPFLAGS '-o ess'
gfortran $FSRCS $CPPFLAGS -o ess || { echo "Compile failed"; exit 1; }

# Model run setting
#OutDir='/media/eweng/HD2/weng/NfixationRuns/'
OutDir='./output/Simulations/'
#.0 .005 .01 .015 .02 .025
# nPFTs = number of PFTs per run; PFTGroup = the PFT combination (PFT 1, 2, 3)
#for nPFTs in 3 2;do
for nPFTs in 2 3;do
if [ "$nPFTs" -eq 3 ]; then
    PFTGroups=("1, 2, 3")
else
    PFTGroups=("1, 2" "1, 3" "2, 3")
fi
for PFTGroup in "${PFTGroups[@]}";do
# R0Nfix only matters when PFT 3 (facultative fixer) is included
if [[ "$PFTGroup" == *3* ]]; then
    R0list=".0 .005 .01 .015 .02 .025"
else
    R0list=".0"
fi
for R0Nfix in $R0list;do
# Run tag without the decimal point, e.g. PFTs123FacuR0N005 (R0Nfix=.005)
groupTag=${PFTGroup//, /}
if [[ "$PFTGroup" == *3* ]]; then
    runTag='PFTs'${groupTag}'FacuR0N'${R0Nfix#.}
else
    runTag='PFTs'${groupTag}
fi

# Directory for the simulation results
DIRECTORY=$OutDir$runTag
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
fp1='./para_files/parameters_BCI_N-Rev1.nml'
echo $fp1

# Run the model with different N deposition
# 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8
for iNdep in 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8;do
   expID="BCI_${runTag}_Ndp${iNdep}"
   echo $expID
   fp2=$DIRECTORY'/parameters_'$expID'.nml'
   sed -e "s#TargetDir#$DIRECTORY#g" \
       -e "s/BCI_NfixTests/$expID/g" \
       -e "s/N_PFTs/$nPFTs/g" \
       -e "s/PFTGroup/$PFTGroup/g" \
       -e "s/R0_Facu/$R0Nfix/g" \
       -e "s/Ndeposition/$iNdep/g" \
       $fp1 > $fp2

   nohup ./ess $fp2 > $expID'.out' 2>&1 &
   echo "$expID submitted."

   sleep 15 # seconds
done
done
done
done

#Clean up compiled files
rm ess
rm *.mod
