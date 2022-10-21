#!/bin/bash

# munge GWAS
ldsc/munge_sumstats.py --sumstats ../frozen_results_files/pw_fetal_sex_gest1.tbl.gz --out pw_fetal_sex_gest_munged --merge-alleles ldscores/w_hm3.snplist --a1 Allele1 --a2 Allele2 --N-col TOTALSAMPLESIZE --snp rsid
ldsc/munge_sumstats.py --sumstats ../frozen_results_files/pw_maternal_sex_gest1.tbl.gz --out pw_maternal_sex_gest_munged --merge-alleles ldscores/w_hm3.snplist --a1 Allele1 --a2 Allele2 --N-col TOTALSAMPLESIZE --snp rsid
ldsc/munge_sumstats.py --sumstats ../frozen_results_files/pw_paternal_sex_gest1.tbl.gz --out pw_paternal_sex_gest_munged --merge-alleles ldscores/w_hm3.snplist --a1 Allele1 --a2 Allele2 --N-col TOTALSAMPLESIZE --snp rsid

# run ldsc
/gpfs/mrc0/projects/Research_Project-MRC158833/rnb203/programs/ldscore/ldsc/ldsc.py --rg pw_fetal_sex_gest_munged.sumstats.gz,pw_maternal_sex_gest_munged.sumstats.gz ----ref-ld-chr ldscores/baselineLD_v1.1/baselineLD. --w-ld-chr ldscores/eur_w_ld_chr/ --out fetal_maternal_sex_gest_ldsc
/gpfs/mrc0/projects/Research_Project-MRC158833/rnb203/programs/ldscore/ldsc/ldsc.py --rg pw_fetal_sex_gest_munged.sumstats.gz,pw_paternal_sex_gest_munged.sumstats.gz ----ref-ld-chr ldscores/baselineLD_v1.1/baselineLD. --w-ld-chr ldscores/eur_w_ld_chr/ --out fetal_paternal_sex_gest_ldsc
/gpfs/mrc0/projects/Research_Project-MRC158833/rnb203/programs/ldscore/ldsc/ldsc.py --rg pw_maternal_sex_gest_munged.sumstats.gz,pw_paternal_sex_gest_munged.sumstats.gz ----ref-ld-chr ldscores/baselineLD_v1.1/baselineLD. --w-ld-chr ldscores/eur_w_ld_chr/ --out maternal_paternal_sex_gest_ldsc
