#!/bin/sh
FSRCS="src/datatypes.F90 \
       src/io_mod.F90 \
       src/soil.F90 \
       src/vegetation.F90 \
       src/BiomeE.F90 \
       src/main.F90"

CPPFLAGS=''
#CPPFLAGS+=' -DScreenOutput'
CPPFLAGS+=' -DDroughtFMT'
CPPFLAGS+=' -DDroughtPaleo'
#CPPFLAGS+=" -DHydro_test"
#CPPFLAGS+=" -DDroughtMu"

echo $FSRCS

#gfortran $FSRCS -o ess -I/opt/local/include -L/opt/local/lib -lnetcdff
#gfortran $FSRCS -o ess -I/Users/eweng/MACPORTS/gcc49-python3/include -L/Users/eweng/MACPORTS/gcc49-python3/lib -lnetcdff
#gfortran src/datatypes.F90 src/io_mod.F90 src/soil.F90 src/vegetation.F90 src/BiomeE.F90 src/main.F90 -DHydro_test -o ess

gfortran $FSRCS $CPPFLAGS -o ess

fparameter='./para_files/parameters_DroughtPaleo.nml'
echo $fparameter
ClimFile='DroughtPaleo_RMA_yrs1514_forcing.csv'
#ClimFile='DroughtPaleo_RMA_forcing.csv'
Run_years='1600'

#for iDraw in {1..1000}; do
for iDraw in {1..10}; do
  runID='Paleo_RMA_'$iDraw
  fp2='./para_files/Paleo_nml/parameters_'$runID'.nml'
  echo "Model run: " $runID
  sed -e "s/Draw_No/$iDraw/g" \
      -e "s/PaleoRunID/$runID/g" \
      -e "s/ForcingFileName/$ClimFile/g" \
      -e "s/RunYears/$Run_years/g" \
      $fparameter > $fp2

  cat $fp2 > ./para_files/input.nml
  ./ess
done
rm ./para_files/input.nml
rm ess
rm esdvm.mod
rm datatypes.mod
rm io_mod.mod
rm soil_mod.mod
rm biomee_mod.mod
