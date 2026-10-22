#!/bin/bash
FSRCS="src/datatypes.F90 \
       src/model_utils.F90 \
       src/io_mod.F90 \
       src/netcdf_io.F90 \
       src/soil.F90 \
       src/vegetation.F90 \
       src/restart_mod.F90 \
       src/BiomeE.F90 \
       src/main.F90"

CPPFLAGS=''
CPPFLAGS+=' -DGlobalRun'
CPPFLAGS+=' -DDO_Climate_VEG'
CPPFLAGS+=' -DUse_GridCellData'
CPPFLAGS+=' -DUse_6hourlyData'
#CPPFLAGS+=' -DZippedNCfiles'
#CPPFLAGS+=' -DHydro_test'
#CPPFLAGS+=' -DSingleTreeTest'
#CPPFLAGS+=' -DScreenOutput'
CPPFLAGS+=' -DZip_outputs'
#CPPFLAGS+=' -DHistCO2'

#NETCDFFLAG='-I/opt/local/include -L/opt/local/lib -lnetcdff'
NETCDFFLAG='-I/usr/local/include -L/usr/local/lib -lnetcdff'

echo $FSRCS
echo $CPPFLAGS
echo $NETCDFFLAG

#gfortran $FSRCS -o ess -I/opt/local/include -L/opt/local/lib -lnetcdff
#gfortran $FSRCS -o ess -I/Users/eweng/MACPORTS/gcc49-python3/include -L/Users/eweng/MACPORTS/gcc49-python3/lib -lnetcdff
#gfortran src/datatypes.F90 src/io_mod.F90 src/soil.F90 src/vegetation.F90 src/BiomeE.F90 src/main.F90 -DHydro_test -o ess
#gfortran -fopenmp $FSRCS $CPPFLAGS -o ess_global -I/opt/local/include -L/opt/local/lib -lnetcdff

gfortran $FSRCS $CPPFLAGS -o ess_global $NETCDFFLAG || { echo "Compile failed"; exit 1; }

# namelist file (Parameter and model setting file)
fparameter='./para_files/parameters_Global_test.nml'

echo $fparameter

# Write to the file that will be read by the model

# Run model
./ess_global $fparameter

rm ess_global
rm *.mod
