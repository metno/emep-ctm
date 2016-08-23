#!/bin/bash

# Minimalistic script for run the Unified EMEP model

# Link the input data
inputdir=.
met=$inputdir/meteo2013
ln -s $inputdir/input/* .   # input files excpet metdata
ln -s $met/DegreeDayFactors.nc .

# Run the model
mpirun $inputdir/code/Unimod

# Clean the links to the input data and remove INPUT.PARA
find . -type l -print0 | xargs -0 rm 

