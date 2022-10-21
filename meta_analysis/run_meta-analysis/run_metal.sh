#!/bin/bash
#SBATCH --export=ALL    # export environment variables (PBS -V)
#SBATCH -D .    # set working directory to . (PBS -d)
#SBATCH --time=48:0:0   # set walltime (PBS -l walltime)
#SBATCH -p mrcq # queue (PBS -q)
#SBATCH -A Research_Project-MRC158833   # research project to run under
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=PW_metal
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# Run meta analyses
/gpfs/mrc0/projects/Research_Project-MRC158833/programs/metal/metal metal_script.txt
gzip -f *.tbl

# Get a list of genome wide significant SNPs
./list_signals.pl --files pw_fetal_sex1.tbl.gz
./list_signals.pl --files pw_fetal_sex_gest1.tbl.gz
./list_signals.pl --files pw_maternal_sex1.tbl.gz
./list_signals.pl --files pw_maternal_sex_gest1.tbl.gz
./list_signals.pl --files pw_paternal_sex1.tbl.gz
./list_signals.pl --files pw_paternal_sex_gest1.tbl.gz

# Generate manhattan and qq plots from meta analyses
./format_for_R.pl --infile pw_fetal_sex1.tbl.gz.filtered.gz --outfile pw_fetal_sex_for_man.gz
./format_for_R.pl --infile pw_fetal_sex_gest1.tbl.gz.filtered.gz --outfile pw_fetal_sex_gest_for_man.gz
./format_for_R.pl --infile pw_maternal_sex1.tbl.gz.filtered.gz --outfile pw_maternal_sex_for_man.gz
./format_for_R.pl --infile pw_maternal_sex_gest1.tbl.gz.filtered.gz --outfile pw_maternal_sex_gest_for_man.gz
./format_for_R.pl --infile pw_paternal_sex1.tbl.gz.filtered.gz --outfile pw_paternal_sex_for_man.gz
./format_for_R.pl --infile pw_paternal_sex_gest1.tbl.gz.filtered.gz --outfile pw_paternal_sex_gest_for_man.gz

module purge
module load R/3.5.1-foss-2018b
Rscript plot_manhattan.R pw_fetal_sex_for_man.gz pw_fetal_sex
Rscript plot_manhattan.R pw_fetal_sex_gest_for_man.gz pw_fetal_sex_gest
Rscript plot_manhattan.R pw_maternal_sex_for_man.gz pw_maternal_sex
Rscript plot_manhattan.R pw_maternal_sex_gest_for_man.gz pw_maternal_sex_gest
Rscript plot_manhattan.R pw_paternal_sex_for_man.gz pw_paternal_sex
Rscript plot_manhattan.R pw_paternal_sex_gest_for_man.gz pw_paternal_sex_gest
