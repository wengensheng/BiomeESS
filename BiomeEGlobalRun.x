#!/bin/sh
FSRCS="src/datatypes.F90 \
       src/io_mod.F90 \
       src/netcdf_io.F90 \
       src/soil.F90 \
       src/vegetation.F90 \
       src/BiomeE.F90 \
       src/BiomeE_global.F90"

CPPFLAGS=''
#CPPFLAGS+="-DHydro_test"
#CPPFLAGS+=' -DSingleTreeTest'
CPPFLAGS+=" -DMergeLowDensityCohorts"
#CPPFLAGS+=" -DScreenOutput"
CPPFLAGS+=" -DGlobalRun"

echo $FSRCS
echo $CPPFLAGS

#gfortran $FSRCS -o ess -I/opt/local/include -L/opt/local/lib -lnetcdff
#gfortran $FSRCS -o ess -I/Users/eweng/MACPORTS/gcc49-python3/include -L/Users/eweng/MACPORTS/gcc49-python3/lib -lnetcdff
#gfortran src/datatypes.F90 src/io_mod.F90 src/soil.F90 src/vegetation.F90 src/BiomeE.F90 src/main.F90 -DHydro_test -o ess

gfortran $FSRCS $CPPFLAGS -o ess_global -I/opt/local/include -L/opt/local/lib -lnetcdff

# namelist file (Parameter and model setting file)
fparameter='./para_files/parameters_Global_test.nml'


echo $fparameter

# Write to the file that will be read by the model
cat $fparameter > ./para_files/input.nml

# Run model
./ess_global

rm ./para_files/input.nml
rm ess_global
rm esdvm.mod
rm datatypes.mod
rm io_mod.mod
rm soil_mod.mod
rm biomee_mod.mod
rm netcdf_io.mod
