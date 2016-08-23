#!/bin/bash

# Minimalistic script for run the Unified EMEP model

# Link the input data
#inputdir=/home/metno/mifasv/OpenSource201310
inputdir=.
ln -s $inputdir/met/*   .   # Driving meteorology
ln -s $inputdir/input/* .   # Other input files

# Define some run parameters
trendyear=2011              # emission year
runlabel1=Base              # short label
runlabel2=Opensource_setup  # long label
startdate="2011 01 01"      # start date (metdata)
  enddate="2011 01 01"      # end date (metdata)

# Put the run parameters in a temporary file
cat > INPUT.PARA << EOF
$trendyear
$runlabel1
$runlabel2
$startdate
$enddate
EOF

# Run the model
mpiexec $inputdir/code/Unimod

# Clean the links to the input data and remove INPUT.PARA
ls $inputdir/met  |xargs rm
ls $inputdir/input|xargs rm
rm INPUT.PARA
