#/bin/sh
gfortran -cpp src/datatypes.F90 src/soil.F90 src/vegetation.F90 src/main.F90 -o ess -I/opt/local/include -L/opt/local/lib -lnetcdff

rm esdvm.mod
rm datatypes.mod
rm soil_mod.mod

./ess

rm ess



