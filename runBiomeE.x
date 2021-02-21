#/bin/sh
#gfortran src/datatypes.F90 src/soil.F90 src/vegetation.F90 src/main.F90 -o ess -I/opt/local/include -L/opt/local/lib -lnetcdff
gfortran src/datatypes.F90 src/soil.F90 src/vegetation.F90 src/main.F90 -o ess -I/Users/eweng/MACPORTS/gcc49-python3/include -L/Users/eweng/MACPORTS/gcc49-python3/lib -lnetcdff
rm esdvm.mod
rm datatypes.mod
rm soil_mod.mod

./ess

rm ess



