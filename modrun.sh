#!/bin/bash
	    
# Minimalistic script for run the EMEP/MSC-W model

# Link the input data
inputdir=.
ln -s $inputdir/input/* .   # Other input files

# Run the model
mpirun $inputdir/code/Unimod

# Clean the links to the input data
find . -type l -delete
