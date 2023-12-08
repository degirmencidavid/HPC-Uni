#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4571mb
#SBATCH --time=00:02:00
#SBATCH --reservation PX425

# NOTE: For OpenMP, we change cpus-per-task and keep ntasks-per-node at 1

# Change into the directory where I submitted the job
cd $SLURM_SUBMIT_DIR 

# Load the modules used when compiling my code
module load GCC

# Set the number of threads to use, this should match
# the number of processors per node (ppn) above.
export OMP_NUM_THREADS=8

# Optional (if using MKL): set the number of thread to use 
# in MKL routines (usually the same as OMP)
# export MKL_NUM_THREADS=$OMP_NUM_THREADS

# Run the executable 5 times with the same settings
./wifi -L 50 -P 0.30 -D 0 -N 5

exit 0
