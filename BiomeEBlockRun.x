#!/bin/bash

# Strings that will be replaced
G2='Grid2'

# namelist file (Parameter and model setting file)
nmlID='Block'$G2
echo "Model run: " $nmlID

# Run model
./ess_global
