#!/bin/bash

#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-15:00 # Runtime in D-HH:MM
#SBATCH -p zorana # Partition to submit to
#SBATCH --mem=100 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o test.out # File to which STDOUT will be written
#SBATCH -e test.err # File to which STDERR will be written
#SBATCH --array=0-19
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=laura.pangx@gmail.com # Email to which notifications will be sent

./a.out 0728_2_4spcs_$SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID 
