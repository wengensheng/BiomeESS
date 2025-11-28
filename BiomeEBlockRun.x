#!/bin/bash

# Strings that will be replaced
G2='Grid2'

# namelist file (Parameter and model setting file)
nmlID='Block'$G2
echo "Model run: " $nmlID

# Write the new nml file (fp2) to the file that will be read by the model

# Run model
./ess_global
