# Code adapted from https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome, November 19th, 2021

import sys
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis


#Read args from command line
args = sys.argv[1:]

plink_prefix_path = args[0]
dnam_bed = args[1]
prefix = args[2]
covariates_file = args[3]
output_dir = args[4]

#Load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(dnam_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

#PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

#Run TensorQTL
cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix, covariates_df=covariates_df, window=500000, output_dir=output_dir)