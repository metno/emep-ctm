#!/bin/bash 

# Minimalistic script to run the Unified EMEP model
#-----------------------------------------------------
# All input variables needed to run the model is now 
# moved into config_emep.nml.   Edit this file for
# iyr_trend, runlabel1, runlabel2, startdate, startdate
# and the link to Meteorology data


# Link the input data
inputdir=.
ln -s $inputdir/input/* .   # Other input files

# Run the model
mpiexec $inputdir/code/Unimod

# Clean the links to the input data. 
ls $inputdir/met  |xargs rm
ls $inputdir/input|xargs rm
