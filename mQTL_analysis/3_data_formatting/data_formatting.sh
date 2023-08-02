#!/bin/sh

ml plink/1.9b_5.2-x86_64
ml StdEnv/2020 r/4.1
ml mugqic/samtools

# Intersect the samples from DNAm ans WGS datasets
mkdir ../formatted_data
awk 'NR==FNR{a[$1];next}($1 in a)' ../Filtered_probes/samples.txt ../final_callset/all_auto.fam > ../formatted_data/intersected_samples.txt

# Format DNAm data
Rscrit --vanilla Format_DNAm_data.R ../Filtered_probes/betas_filter_confidence_variance.RData ../formatted_data/intersected_samples.txt ../formatted_data/format_DNAm_data.txt ../formatted_data/betas.bed
sort -k 1,1 -k2,2n ../formatted_data/betas.bed > ../formatted_data/betas_sorted.bed
bgzip ../formatted_data/betas_sorted.bed && tabix -p bed ../formatted_data/betas_sorted.bed.gz

# Format WGS data
plink --bfile ../final_callset/all_auto --make-bed --keep-fam ../formatted_data/intersected_samples.txt --geno 0.2 --mac 10 --hwe 0.000001 --out ../formatted_data/all_auto

# Format covariates
Rscript --vanilla Format_covariates.R ../formatted_data/intersected_samples.txt ../PACE_pre_processing/estimated_cell_types.tsv ../pheno_sans_dup_avec_predg_2021_07_07.csv ../Relatedness/geneticPCs.tsv  ../formatted_data/Format_covariates.txt ../formatted_data/covariates.bed