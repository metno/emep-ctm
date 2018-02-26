#!/bin/bash

# Minimalistic script for run the Unified EMEP model
GRID=EECCA
NLEV=20lev
cd ~/work/EMEP_MSC-W_model.rv4.17.OpenSource/Base_${GRID}_${NLEV}

# Run the model
mpiexec ../code/Unimod
