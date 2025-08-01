#!/bin/sh
FSRCS="src/datatypes.F90 \
       src/io_mod.F90 \
       src/soil.F90 \
       src/vegetation.F90 \
       src/BiomeE.F90 \
       src/main.F90"

CPPFLAGS=''
#CPPFLAGS+="-DHydro_test"
#CPPFLAGS+=' -DSingleTreeTest'
CPPFLAGS+=" -DMergeLowDensityCohorts"
CPPFLAGS+=" -DScreenOutput"

echo $FSRCS

#gfortran $FSRCS -o ess -I/opt/local/include -L/opt/local/lib -lnetcdff
#gfortran $FSRCS -o ess -I/Users/eweng/MACPORTS/gcc49-python3/include -L/Users/eweng/MACPORTS/gcc49-python3/lib -lnetcdff
#gfortran src/datatypes.F90 src/io_mod.F90 src/soil.F90 src/vegetation.F90 src/BiomeE.F90 src/main.F90 -DHydro_test -o ess

gfortran $FSRCS $CPPFLAGS -o ess

# namelist file (Parameter and model setting file)
#fparameter='./para_files/parameters_BCI_EG_DC.nml'
fparameter='./para_files/parameters_ORNL_test.nml'


echo $fparameter

# Write to the file that will be read by the model
cat $fparameter > ./para_files/input.nml

# Run model
./ess

rm ./para_files/input.nml
rm ess
rm esdvm.mod
rm datatypes.mod
rm io_mod.mod
rm soil_mod.mod
rm biomee_mod.mod
