#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4571mb
#SBATCH --time=00:02:00
#SBATCH --reservation PX425

# NOTE: For MPI, we change ntasks-per-node and keep  cpus-per-task at 1

# Change into the directory where I submitted the job
cd $SLURM_SUBMIT_DIR 

# Load the modules used when compiling my code
module load GCC

# Note we use srun rather than mpirun here
srun ./wifi -L 50 -P 0.30 -D 0 -N 5

exit 0
