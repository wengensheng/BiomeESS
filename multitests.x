#/bin/sh
gfortran src/datatypes.F90 src/soil.F90 src/vegetation.F90 src/main_multi_tests.F90 -o ess2
./ess2

rm ess2
rm esdvm.mod
rm datatypes.mod
rm soil_mod.mod

