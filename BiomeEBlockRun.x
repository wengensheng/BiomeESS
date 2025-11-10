#!/bin/bash

# Strings that will be replaced
G1='Grid1'
G2='Grid2'

# namelist file (Parameter and model setting file)
nmlID='Block'$G2
echo "Model run: " $nmlID
fp1='./para_files/parameters_BlockRun.nml'
fp2='./para_files/parameters_'$nmlID'.nml'
sed -e "s/StartGrid/$G1/g" \
    -e "s/EndGrid/$G2/g" \
    $fp1 > $fp2

# Write the new nml file (fp2) to the file that will be read by the model
cat $fp2 > ./para_files/input.nml

# Run model
#nohup ./ess_global > 'GridRun_'$G2'.out' 2>&1 &
./ess_global
