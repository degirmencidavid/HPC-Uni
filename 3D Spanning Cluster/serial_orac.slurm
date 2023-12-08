#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4571mb
#SBATCH --time=00:02:00
#SBATCH --reservation PX425

# Change into the directory where I submitted the job
cd $SLURM_SUBMIT_DIR 

export OMP_NUM_THREADS=1
./wifi -L 50 -P 0.30 -D 0 -N 5

exit 0
