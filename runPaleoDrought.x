#!/bin/bash
set -euo pipefail

# --------------------------------------------------
# Source files
# --------------------------------------------------
FSRCS="src/datatypes.F90 \
       src/io_mod.F90 \
       src/soil.F90 \
       src/vegetation.F90 \
       src/BiomeE.F90 \
       src/main.F90"

CPPFLAGS=''
CPPFLAGS+=' -DDroughtFMT'
CPPFLAGS+=' -DDroughtPaleo'

# --------------------------------------------------
# Compile ONCE
# --------------------------------------------------
echo "Compiling BiomeE..."
gfortran $FSRCS $CPPFLAGS -o BiomeE

# --------------------------------------------------
# Site / PFT mapping
# --------------------------------------------------
Sites=("RMA" "SJC" "DCK")
PFT_list=("7, 0" "4, 1" "0, 1")

Run_years='1699'
fp_template='./para_files/parameters_DroughtPaleo.nml'

BaseDir='/media/eweng/HD2/weng/PaleoTests'

# --------------------------------------------------
# Loop over sites
# --------------------------------------------------
for idx in "${!Sites[@]}"; do

  Site=${Sites[$idx]}
  PFTs="${PFT_list[$idx]}"
  SiteDir="${BaseDir}/Paleo${Site}"

  echo "===================================="
  echo "Running Site: $Site"
  echo "PFTs: $PFTs"
  echo "Output directory: $SiteDir"
  echo "===================================="

  # ------------------------------------------------
  # Check if directory exists
  # ------------------------------------------------
  if [ ! -d "$SiteDir" ]; then
      echo "Directory $SiteDir does not exist. Creating it..."
      mkdir -p "$SiteDir"
      if [ $? -ne 0 ]; then
          echo "ERROR: Failed to create directory $SiteDir"
          exit 1
      fi
  else
      echo "Directory $SiteDir already exists."
  fi

  # ------------------------------------------------
  # Loop over stochastic draws
  # ------------------------------------------------
  for iDraw in {1..5}; do

    runID="Paleo_${Site}_${iDraw}"
    paramFile="${SiteDir}/parameters_${runID}.nml"

    echo "Running: $runID"

    sed -e "s/Draw_No/${iDraw}/g" \
        -e "s/PaleoRunID/${runID}/g" \
        -e "s/SiteID/${Site}/g" \
        -e "s/SiteSP/${PFTs}/g" \
        -e "s/RunYears/${Run_years}/g" \
        -e "s#TargetDir#${SiteDir}#g" \
        "$fp_template" > "$paramFile"

    cp "$paramFile" ./para_files/input.nml
    ./BiomeE
    rm -f ./para_files/input.nml

  done
done

# --------------------------------------------------
# Cleanup
# --------------------------------------------------
rm -f *.mod
echo "All site simulations completed successfully."
