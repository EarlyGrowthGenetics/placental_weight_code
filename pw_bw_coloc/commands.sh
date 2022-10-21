#!/bin/bash

#TODO convert b37 to b38
#TODO lift over (I'd done b38 positions but lifted from b37 to 38!!!)
#TODO rewrite extract_and_join.pl to lift the right way

# extract positions from deCODE SNPs to lift over to b37
total=`wc -l decode_loci_positions_b37 | awk '{print $1}'`
dfetal="~/decode_BW_gwas/Birthweight2021.gz"
dmaternal="~/decode_BW_gwas/Birthweight_offspring_mothers2021.gz"
pfetal="../frozen_results_files/pw_fetal_sex_gest_1.gz"
pmaternal="../frozen_results_files/pw_maternal_sex_gest_1.gz"
ppaternal="../frozen_results_files/pw_paternal_sex_gest.gz"
for i in $(seq 1 $total)
do
	# get the locus information
	IFS=$'\n'
	line1=`head -n $i decode_loci_positions_b37 | tail -n 1`
	line2=`head -n $i decode_loci_positions_b38 | tail -n 1`
	IFS=$'\t' read -r -a pos37 <<< $line1
	IFS=$'\t' read -r -a pos38 <<< $line2

	# decide which decode file we need
	if [[ ${pos37[4]} -eq "Fetal" ]]
	then
		decode=$pfetal
	elif [[ ${pos37[4]} -eq "Maternal" ]]
	then
		decode=$pmaternal
	else
		decode=$dpaternal
	fi

	# extract the locus information from decode files
	zcat $decode | awk -v chr=${pos37[0]} -v pos=${pos37[1]} '$17==chr && $18>pos-500000 && $18<pos+500000{print "chr"$17"\t"$18"\t"$18}' > locus_$i
	liftOver locus_${i} hg19ToHg38.over.chain.gz locus_${i}_lifted locus_${i}_failed
done

# join PW and deCODE summary stats for SNP of interest
./extract_and_join.pl

# run coloc
echo "" > results_joined
for i in $(seq 1 $total)
do
	Rscript run_coloc.R $i
	tail -n 1 locus_${i}_out >> results_joined
done
