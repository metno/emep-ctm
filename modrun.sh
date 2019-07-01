#!/bin/bash

### job options for Slurm/sbatch
#SBATCH --job-name=emepctm
#SBATCH --output=%x.out --error=%x.out
#SBATCH --nodes=2 --ntasks-per-node=32 --time=5:00:00

### Minimalistic script for run the Unified EMEP model

# working directory
cd ~/work/EMEP_MSC-W_model.rv4.33.OpenSource/Base

# run the model
#mpiexec ../code/emepctm # or
#mpirun  ../code/emepctm # or
mpprun   ../code/emepctm # depending on the HPC/queue



