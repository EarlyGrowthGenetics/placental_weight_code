#!/bin/bash
#SBATCH --export=ALL	# export environment variables (PBS -V)
#SBATCH -D .	# set working directory to . (PBS -d)
#SBATCH --time=48:0:0	# set walltime (PBS -l walltime)
#SBATCH -p mrcq	# queue (PBS -q)
#SBATCH -A Research_Project-MRC158833	# research project to run under
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-4

i=$SLURM_ARRAY_TASK_ID 

./merge.pl --locus $i --analysis gest
module purge
module load R/4.0.0-foss-2020a
Rscript run_coloc.R $i gest
