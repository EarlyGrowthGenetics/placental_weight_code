#!/bin/sh

ml mugqic/bcftools
ml plink/2.00a2.3_x86_64
ml StdEnv/2020 r/4.1

# Keep only SNPs and split multiallelic sites
mkdir ../snps_split_multiallelic
for chr in {1..22}; do
    bcftools norm -m -any -Ou children_clean_$chr.vcf.gz | bcftools view -v snps -Oz -o ../snps_split_multiallelic/wgs_$chr.vcf.gz
done

# Filter depth and callrate
mkdir ../depth_callrate
for chr in {1..22}; do
    bcftools filter -e "FORMAT/DP<10" -S . -Ou ../snps_split_multiallelic/wgs_$chr.vcf.gz | bcftools filter -e "F_MISSING>0.2" -Oz -o ../depth_callrate/wgs_$chr.vcf.gz
done

# Convert to plink format and concatenate
mkdir ../plink
for chr in {1..22}; do
    plink2 --vcf ../depth_callrate/wgs_$chr.vcf.gz --make-bed --snps-only --autosome --set-missing-var-ids @:#[b37]\$r,\$a --new-id-max-allele-len 60 --memory 5000 --out ../plink/wgs_$chr
done

## Create list of plink files 
for file in $(ls -d ../plink/*); do
    echo "${{file%.*}}" >> ../plink/files_list.txt;
done
sort ../plink/files_list.txt | uniq > ../plink/files_list_sorted.txt


ml plink/1.9b_5.2-x86_64
## Merge plink files
plink --merge-list ../plink/files_list_sorted.txt --make-bed --memory 5000 --out ../plink/all_auto --update-sex ../sexInfoFormatted.tab

# Missingness and heterozygosity
mkdir ../missing_heterozygosity

## Missingness
plink --bfile ../plink/all_auto --missing --out ../missing_heterozygosity/all_auto
sed -i 's/[ ]\+/\t/g' ../missing_heterozygosity/all_auto.imiss
sed -i 's/^\t//' ../missing_heterozygosity/all_auto.imiss

## Heterozygocity
plink --bfile ../plink/all_auto --het --out ../missing_heterozygosity/all_auto
sed -i 's/[ ]\+/\t/g' ../missing_heterozygosity/all_auto.het
sed -i 's/^\t//' ../missing_heterozygosity/all_auto.het

## Plot missingness and heterozygosity
Rscript --vanilla ../missing_heterozygosity/all_auto.imiss ../missing_heterozygosity/all_auto.het ../missing_heterozygosity/miss.png ../missing_heterozygosity/het.png
### Manually check if any heterozygosity or missingness outliers and report them in ../missing_heterozygosity/outliers.txt

# Make stringent callset for population analyzes
mkdir ../stringent_callset

plink --bfile ../plink/all_auto --make-bed --geno 0.01 --maf 0.05 --hwe 0.0001 --remove-fam ../missing_heterozygosity/outliers.txt --out ../stringent_callset/all_auto
plink --bfile ../stringent_callset/all_auto --recode vcf-iid bgz --memory 30000 --threads 10 --out ../stringent_callset/all_auto_filtered
bcftools index -t ../stringent_callset/all_auto_filtered.vcf.gz

# Merge with 1000 Genomes project data
mkdir ../1000G_merged
1kg_vcf=../1kg_phase1_all_SNPsOnly_maf05_cr99.vcf.gz

bcftools isec -n=2 -p ../1000G_merged -Oz ../stringent_callset/all_auto_filtered.vcf.gz $1kg_vcf
bcftools merge -m none -Oz -o ../1000G_merged/Gen3G_stringent_1000G.vcf.gz ../1000G_merged/0000.vcf.gz ../1000G_merged/0001.vcf.gz
bcftools index --threads 10 ../1000G_merged/Gen3G_stringent_1000G.vcf.gz

# Convert merged vcf to GDS
Rscript --vanilla convertGDS.R ../1000G_merged/Gen3G_stringent_1000G.vcf.gz ../1000G_merged/Gen3G_stringent_1000G.gds

# Estimate ancestry in Gen3G based on 1000G data with the GENESIS suite
mkdir ../Ancestry

Rscript --vanilla Estimate_ancestry.R ../Ancestry ../1000G_merged/Gen3G_stringent_1000G.gds ../Ancestry/Estimate_ancestry.log ../igsr-1000_genomes_phase_1_release.csv ../Ancestry/Outliers/outliers.txt
### Manually validate ancestry outliers and report those to remove in ../Ancestry/outliers_to_remove.txt

# Estimate relatedness in Gen3G with the GENESIS suite
mkdir ../Relatedness

Rscript --vanilla Estimate_relatedness.R ../Relatedness ../1000G_merged/Gen3G_stringent_1000G.gds ../Ancestry/outliers_to_remove.txt ../Relatedness/estimate_relatedness.log ../Relatedness/geneticPCs.tsv ../Relatedness/related.tsv
### Manually check the related pairs, select which one to remove based on sewuencing QC and report them in ../Relatedness/related_to_remove.txt

# Make final callset
mkdir ../final_callset

cat ../Ancestry/outliers_to_remove.txt ../Relatedness/related_to_remove.txt > ../final_callset/to_remove.txt
grep "*" ../plink/all_auto.bim | cut -f2 > ../final_callset/starGeno
plink --bfile ../plink/all_auto --make-bed --remove ../final_callset/to_remove.txt --exclude ../final_callset/starGeno --geno 0.2 --mac 10 --hwe 0.000001 --out ../final_callset/all_auto