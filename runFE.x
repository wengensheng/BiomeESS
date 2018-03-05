#/bin/sh
gfortran src/datatypes.F90 src/soil.F90 src/vegetation.F90 src/main.F90 -o ess
./ess

rm ess
rm esdvm.mod
rm datatypes.mod
rm soil_mod.mod

