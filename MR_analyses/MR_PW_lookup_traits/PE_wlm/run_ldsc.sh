#!/bin/bash

/slade/home/rnb203/programs/ldscore/ldsc/munge_sumstats.py --sumstats fetal_preeclampsia_filtered.gz --out fetal_PE_munged --merge-alleles /slade/home/rnb203/programs/ldscore/ldscores/w_hm3.snplist --a1 Allele1 --a2 Allele2 --N-cas 4630 --N-con 373345 --snp MarkerName
/slade/home/rnb203/programs/ldscore/ldsc/munge_sumstats.py --sumstats maternal_preeclampsia_filtered.gz --out maternal_PE_munged --merge-alleles /slade/home/rnb203/programs/ldscore/ldscores/w_hm3.snplist --a1 Allele1 --a2 Allele2 --N-cas 7219 --N-con 155660 --snp MarkerName
/slade/home/rnb203/programs/ldscore/ldsc/ldsc.py --h2 fetal_PE_munged.sumstats.gz,maternal_PE_munged.sumstats.gz ----ref-ld-chr /slade/home/rnb203/programs/ldscore/ldscores/baselineLD_v1.1/baselineLD. --w-ld-chr /slade/home/rnb203/programs/ldscore/ldscores/eur_w_ld_chr/ --out PE_overlap
