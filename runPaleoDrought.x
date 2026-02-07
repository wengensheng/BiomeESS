#!/bin/bash
FSRCS="src/datatypes.F90 \
       src/io_mod.F90 \
       src/soil.F90 \
       src/vegetation.F90 \
       src/BiomeE.F90 \
       src/main.F90"

CPPFLAGS=''
CPPFLAGS+=' -DScreenOutput'
#CPPFLAGS+=' -DDO_Climate_VEG'
CPPFLAGS+=' -DDroughtFMT'
CPPFLAGS+=' -DDroughtPaleo'
#CPPFLAGS+=' -DDBEN_runs'
#CPPFLAGS+=" -DHydro_test"
#CPPFLAGS+=" -DDroughtMu"

echo $FSRCS

#gfortran $FSRCS -o ess -I/opt/local/include -L/opt/local/lib -lnetcdff
#gfortran $FSRCS -o ess -I/Users/eweng/MACPORTS/gcc49-python3/include -L/Users/eweng/MACPORTS/gcc49-python3/lib -lnetcdff
#gfortran src/datatypes.F90 src/io_mod.F90 src/soil.F90 src/vegetation.F90 src/BiomeE.F90 src/main.F90 -DHydro_test -o ess
#gfortran $FSRCS $CPPFLAGS -o BiomeE_paleo -I/usr/local/include -L/usr/local/lib -lnetcdff
gfortran $FSRCS $CPPFLAGS -o BiomeE_paleo

runTag='PaleoDrought'
#DIRECTORY='./output/'$runTag
DIRECTORY='/media/eweng/HD2/weng/'$runTag

# Check if the directory exists. If not, create it.
if [ ! -d "$DIRECTORY" ]; then
    echo "Directory $DIRECTORY does not exist. Creating it now..."
    mkdir -p "$DIRECTORY"
    if [ $? -eq 0 ]; then
        echo "Directory $DIRECTORY created successfully."
    else
        echo "Failed to create directory $DIRECTORY."
        exit 1 # Exit with an error code if creation fails
    fi
else
    echo "Directory $DIRECTORY already exists."
fi

fp1='./para_files/parameters_DroughtPaleo.nml'
echo $fp1
ClimFile='DroughtPaleo_RMA_yrs1514_forcing.csv'
#ClimFile='DroughtPaleo_RMA_forcing.csv'
Run_years='1699'

#for iDraw in {1..200}; do
for iDraw in {1..1}; do
  runID='Paleo_RMA_'$iDraw
  fp2=$DIRECTORY'/parameters_'$runID'.nml'
  echo "Model run: " $runID
  sed -e "s/Draw_No/$iDraw/g" \
      -e "s/PaleoRunID/$runID/g" \
      -e "s/ForcingFileName/$ClimFile/g" \
      -e "s/RunYears/$Run_years/g" \
      -e "s#TargetDir#$DIRECTORY#g" \
      $fp1 > $fp2

  cat $fp2 > ./para_files/input.nml
  ./BiomeE_paleo
done

rm ./para_files/input.nml
rm BiomeE_paleo
rm *.mod
