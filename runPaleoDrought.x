#!/bin/bash
set -euo pipefail

FSRCS="src/datatypes.F90 \
       src/io_mod.F90 \
       src/soil.F90 \
       src/vegetation.F90 \
       src/BiomeE.F90 \
       src/main.F90"

CPPFLAGS=''
CPPFLAGS+=' -DDroughtFMT'
CPPFLAGS+=' -DDroughtPaleo'

echo "Compiling BiomeE..."
gfortran $FSRCS $CPPFLAGS -o BiomeE
EXE="$(readlink -f ./BiomeE)"   # absolute path to executable

Sites=("SJC" "RMA" "DCK")
PFT_list=("4, 1" "7, 0" "0, 1")

Run_years='1699'
fp_template='./para_files/parameters_DroughtPaleo.nml'
BaseDir='/media/eweng/HD2/weng/PaleoTests'

# -----------------------------
# User-controlled draw range
# -----------------------------
DRAW_FIRST=1
DRAW_LAST=200

# -----------------------------
# Batch / concurrency settings
# -----------------------------
BATCH_SIZE=20   # "20 draws a time"
MAX_JOBS=20     # max concurrent runs (<= 32 cores is fine)

# --- sanity checks ---
if [ "$DRAW_FIRST" -lt 1 ]; then
  echo "ERROR: DRAW_FIRST must be >= 1"
  exit 1
fi
if [ "$DRAW_LAST" -lt "$DRAW_FIRST" ]; then
  echo "ERROR: DRAW_LAST must be >= DRAW_FIRST"
  exit 1
fi

wait_for_slot () {
  while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do
    sleep 1
  done
}

for idx in "${!Sites[@]}"; do
  Site=${Sites[$idx]}
  PFTs="${PFT_list[$idx]}"

  SiteDir="${BaseDir}/Paleo${Site}"

  echo "===================================="
  echo "Running Site: $Site | PFTs: $PFTs"
  echo "Draws: $DRAW_FIRST to $DRAW_LAST"
  echo "Output directory: $SiteDir"
  echo "===================================="

  # ---- keep your directory exists check ----
  if [ ! -d "$SiteDir" ]; then
      echo "Directory $SiteDir does not exist. Creating it..."
      mkdir -p "$SiteDir"
      if [ $? -ne 0 ]; then
          echo "Failed to create directory $SiteDir."
          exit 1
      fi
  else
      echo "Directory $SiteDir already exists."
  fi

  # ---- batches of 20 within [DRAW_FIRST, DRAW_LAST] ----
  for ((batch_start=DRAW_FIRST; batch_start<=DRAW_LAST; batch_start+=BATCH_SIZE)); do
    batch_end=$((batch_start + BATCH_SIZE - 1))
    if [ "$batch_end" -gt "$DRAW_LAST" ]; then
      batch_end=$DRAW_LAST
    fi

    echo "------------------------------------"
    echo "Launching draws $batch_start to $batch_end (up to $MAX_JOBS concurrent)"
    echo "------------------------------------"

    # launch this batch
    for ((iDraw=batch_start; iDraw<=batch_end; iDraw++)); do
      runID="Paleo_${Site}_${iDraw}"
      paramFile="${SiteDir}/parameters_${runID}.nml"
      logFile="${SiteDir}/log_${runID}.txt"

      sed -e "s/Draw_No/${iDraw}/g" \
          -e "s/PaleoRunID/${runID}/g" \
          -e "s/SiteID/${Site}/g" \
          -e "s/SiteSP/${PFTs}/g" \
          -e "s/RunYears/${Run_years}/g" \
          -e "s#TargetDir#${SiteDir}#g" \
          "$fp_template" > "$paramFile"

      (
        "$EXE" "$paramFile" > "$logFile" 2>&1
      ) &

      wait_for_slot
    done

    # wait for all runs in this batch to finish
    wait
    echo "Batch $batch_start-$batch_end finished."
    echo ""
  done
done

rm -f *.mod
echo "All simulations completed."
