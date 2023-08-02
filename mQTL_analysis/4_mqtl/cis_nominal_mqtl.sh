#!/bin/sh

python/3.11.2

# Run tensorQTL using python. One might need to first seperate the data for each chromosome, and even chromosomes into chuncks for the large ones.
mkdir ../mqtl
python3 cis_nominal_mqtl.py ../formatted_data/all_auto ../formatted_data/betas_sorted.bed.gz "2021_07_cis_nominal" ../formatted_data/covariates.bed ../mqtl