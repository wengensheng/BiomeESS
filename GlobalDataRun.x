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

Lon1=(1   121 181 241 361 421 481 541 601)
Lon2=(120 180 240 360 420 480 540 600 720)
# namelist file (Parameter and model setting file)
fp1='./para_files/parameters_GlobalData.nml'
echo $fp1

for iB in "${!Lon2[@]}"; do
  if [ "${Lon1[$iB]}" -lt "${Lon2[$iB]}" ]; then
    runID='Lon_'${Lon2[$iB]}
    fp2='./para_files/parameters_'$runID'.nml'

    echo "Block ${Lon1[$iB]}-${Lon2[$iB]}"
    echo $fp2
    echo "Model run: " $runID
    sed -e "s/LonStart/${Lon1[$iB]}/g" \
        -e "s/LonEnd/${Lon2[$iB]}/g" \
        -e "s/GlobalVegGridList/VegList$runID/g" \
        $fp1 > $fp2

    echo "Run Longitude ${Lon1[$iB]}-${Lon2[$iB]}"
    cat $fp2 > ./para_files/input.nml

    # Run model
    ./ess_global

    rm ./para_files/input.nml
  fi
done


#rm ess_global
rm *.mod
