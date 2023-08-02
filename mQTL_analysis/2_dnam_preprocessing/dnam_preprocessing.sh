#!/bin/sh

ml StdEnv/2020 r/4.1

# Pre-processing of DNAm data based on the PACEAnalysis pipeline
mkdir ../PACE_pre_processing
Rscript --vanilla pace_preprocessing.R ../RGset_ext_combined.RData ../PACE_pre_processing ../pheno_sans_dup_avec_predg_2021_07_07.csv

# Filtering of the probes based on cross-reactivity and SNPs
mkdir ../Filtered_probes
Rscript --vanilla Filter_cpgs_confidence.R ../PACE_pre_processing/Gen3G_20210709_PreprocessedBetas_nooutliers.RData ../13059_2016_1066_MOESM1_ESM.csv ../Filtered_probes/log_confidence.txt ../Filtered_probes/betas_filter_confidence.RData

# Filtering of the probes based on variance
Rscript --vanilla Filter_cpgs_low_variance.R ../Filtered_probes/betas_filter_confidence.RData 0.00001 ../Filtered_probes/log_variance.txt ../Filtered_probes/betas_filter_confidence_variance.RData ../Filtered_probes/samples.txt ../Filtered_probes/variance_distributions.pdf