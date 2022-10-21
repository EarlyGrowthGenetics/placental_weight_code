#!/bin/bash

# extract SNPs required from PW results
# Traits are:
# ins_sec
# ins_res
# Height
# FG
# SBP
# DBP

# generate list of all SNPs we need
printf "ins_sec\nins_res\nHeight\nFG\nSBP\nDBP" > trait_names
grep -wf trait_names snp_info.txt | cut -f 2,3,4 | sort -k1,1 -k2,2 | uniq > SNP_rsids_all_traits	# pull out only unique SNPs, but this assumes alleles are listed the same way around for all traits so there may be duplicate SNPs in this list
# set header for extracted file
printf "Trait\tSNP\tA1\tBetaYG\tseBetaYG\tP_YG\n" > SNPs_all_traits_extracted
# extract SNPs from meta-analysis
./extract_snps.pl --infile ..//frozen_meta_analysis_files/pw_fetal_sex.gz --outfile SNPs_all_traits_extracted --snps SNP_rsids_all_traits --trait PW_fetal_sex
./extract_snps.pl --infile ../frozen_meta_analysis_files/pw_maternal_sex.gz --outfile SNPs_all_traits_extracted --snps SNP_rsids_all_traits --trait PW_maternal_sex
./extract_snps.pl --infile ../frozen_meta_analysis_files/pw_paternal_sex.gz --outfile SNPs_all_traits_extracted --snps SNP_rsids_all_traits --trait PW_paternal_sex
./extract_snps.pl --infile ../frozen_meta_analysis_files/pw_fetal_sex_gest.gz --outfile SNPs_all_traits_extracted --snps SNP_rsids_all_traits --trait PW_fetal_sex_gest
./extract_snps.pl --infile ../frozen_meta_analysis_files/pw_maternal_sex_gest.gz --outfile SNPs_all_traits_extracted --snps SNP_rsids_all_traits --trait PW_maternal_sex_gest
./extract_snps.pl --infile ../frozen_meta_analysis_files/pw_paternal_sex_gest.gz --outfile SNPs_all_traits_extracted --snps SNP_rsids_all_traits --trait PW_paternal_sex_gest
# extract SNPs from WLM
./extract_snps_wlm.pl --infile ../frozen_meta_analysis_files/wlm_sex.gz --outfile SNPs_all_traits_extracted --snps SNP_rsids_all_traits --trait PW_fetal_sex_wlm --beta 18 --se 19
./extract_snps_wlm.pl --infile ../frozen_meta_analysis_files/wlm_sex.gz --outfile SNPs_all_traits_extracted --snps SNP_rsids_all_traits --trait PW_maternal_sex_wlm --beta 20 --se 21
./extract_snps_wlm.pl --infile ../frozen_meta_analysis_files/wlm_sex.gz --outfile SNPs_all_traits_extracted --snps SNP_rsids_all_traits --trait PW_paternal_sex_wlm --beta 22 --se 23
./extract_snps_wlm.pl --infile ../frozen_meta_analysis_files/wlm_sex_gest.gz --outfile SNPs_all_traits_extracted --snps SNP_rsids_all_traits --trait PW_fetal_sex_gest_wlm --beta 18 --se 19
./extract_snps_wlm.pl --infile ../frozen_meta_analysis_files/wlm_sex_gest.gz --outfile SNPs_all_traits_extracted --snps SNP_rsids_all_traits --trait PW_maternal_sex_gest_wlm --beta 20 --se 21
./extract_snps_wlm.pl --infile ../frozen_meta_analysis_files/wlm_sex_gest.gz --outfile SNPs_all_traits_extracted --snps SNP_rsids_all_traits --trait PW_paternal_sex_gest_wlm --beta 22 --se 23

# run the MR script for both
module purge
module load R
for name in `cat trait_names`
do
  mkdir -p $name
  Rscript general_IV_script_plots.R $name snp_info.txt SNPs_all_traits_extracted $name
done

rm trait_names
