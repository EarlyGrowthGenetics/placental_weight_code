#!/bin/bash

cd PE_wlm/
zcat fetal_preeclampsia_filtered.gz | cut -f 1-3,8-10 | gzip -c > fetal_preeclampsia_filtered_tmp.gz
zcat maternal_preeclampsia_filtered.gz | cut -f 1-3,8-10 | gzip -c > maternal_preeclampsia_filtered_tmp.gz
./parent_offspring_wlm/apply_wlm.pl --child_file fetal_preeclampsia_filtered_tmp.gz --mother_file maternal_preeclampsia_filtered_tmp.gz --out_file PE_wlm.gz --cm 0.4445
rm fetal_preeclampsia_filtered_tmp.gz maternal_preeclampsia_filtered_tmp.gz
loadR
Rscript apply_wlm.R
cd ../

./extract_snps.pl --infile PE_wlm/maternal_preeclampsia_filtered.gz --outfile SNPs_all_traits_extracted --snps SNP_rsids_all_traits --trait PE_fetal --rsid 0 --ea 1 --oa 2 --beta 7 --se 8 --p 9
./extract_snps_PE_liability.pl --infile PE_wlm/PE_wlm_adjusted.txt.gz --outfile SNPs_all_traits_extracted --snps SNP_rsids_all_traits --trait PE_fetal_wlm_liability --rsid 0 --ea 1 --oa 2 --beta 9 --se 10 --p 11

for name in PW_gest_fetal PW_gest_fetal_wlm_trio PW_gest_fetal_wlm_pair
do
	Rscript ../general_IV_script_plots.R $name ../snp_info.txt SNPs_all_traits_extracted $name
done
