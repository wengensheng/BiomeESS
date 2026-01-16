#!/bin/bash

# Tag a block ID to the model run
G2='Grid2'
nmlID='Block'$G2
echo "Model run: " $nmlID

# Run model
./ess_global
