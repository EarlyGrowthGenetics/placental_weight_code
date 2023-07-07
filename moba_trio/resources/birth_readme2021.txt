There are five files with GWAS results, three based on a meta-analysis of Iceland and EGG/UKB data,

Birthlength.txt
Birthweigh.txt
Birthweigh_offspring_mothers.txt

and two based on Icelandic data only

Birthweigh_offspring_fathers.txt
Poinderal_index.txt

All the files includ the following columns

Chr		Chromosome
Pos		Position in Human Build 38
rsID		rs-ID if known
A0		Reference allele
A1		Effect allele
IS-frq	Frequency of effect allele in the Icelandic data
IS-info	Imputation information for the variant in the Icelandic data
Beta-A1	Effect estimate in the Icelandic analysis or in the meta-anlysis of Iceland and EGG/UKB
P		P-value from the the Icelandic analysis or in the meta-anlysis of Iceland and EGG/UKB

Additional columns for the three files with meta-analysis results

EGG-frq	Frequency of the effect allele in the EGG/UKB data
I2		I2 hetereogeneity statistics from the meta-analysis
P-het		Heterogeneitiy P value from the meta-analysis
