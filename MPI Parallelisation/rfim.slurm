#!/bin/bash
# Orac submission script for PX425 2021 Assignment 4
# Runs on the "px425" queue which is a 2-node reservation

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3882mb
#SBATCH --time=00:01:00
#SBATCH --reservation PX425

# Some comments on resource requests
# ==================================
#
# There are two ways to request resources on orac:
# One can either use nodes=q ntasks-per-node=p and be 
# allocated p cores on q physical nodes, or one
# can just ask for ntasks=X, and get X processors
# spread over any physical nodes. However, in this
# assignment, we have a reservation of 2 nodes, so
# as long as you use X<=28 they will often all be on
# one node, but it is possible they will be spread out
# over both of the nodes.
# 
# NB: if the system gets busy, you may experience
# some variability of timings larger X since your
# job might be split over both nodes. If that appears
# to be the case, and you need accurate timings, try
# the nodes=1:ntasks-per-node=X option instead.

# Change into the directory where I submitted the job
cd $SLURM_SUBMIT_DIR 

# Load the modules used when compiling my code
module load intel
module load impi

# Note we use srun rather than mpirun here
srun ./rfim 

exit 0
