#!/bin/bash

for i in {1..23}
do
	./generate_files.pl --analysis sex --outfile sex/locus_${i}.gz --snps loci_sex --snp_num ${i}
	gunzip sex/locus_${i}.gz
	chr=`head -n $i loci_sex | tail -n 1 | cut -f 1`
	pos=`head -n $i loci_sex | tail -n 1 | cut -f 2`
	from=$(( $pos-1000000 ))
	to=$(( $pos+1000000 ))
	echo $chr $from $to "set1" > sex/range_${i}
	/gpfs/mrc0/projects/Research_Project-MRC158833/programs/plink2/plink2 --bfile /gpfs/mrc0/projects/Research_Project-MRC158833/UKBiobank/500K_Genetic_data/imputed_data_best_guess/ukb_imp_chr${chr}_v3_best_guess_WB_340K --chr $chr --extract range sex/range_${i} --make-bed --out sex/plink_UKBB_snp_${i}
	/gpfs/mrc0/projects/Research_Project-MRC158833/programs/gcta/gcta64 --bfile sex/plink_UKBB_snp_${i} --chr ${chr} --maf 0.001 --cojo-file sex/locus_${i} --cojo-slct --out sex/cojo_chr${i}_340k --thread-num 16
	gzip sex/locus_${i}
done
