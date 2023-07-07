
##
#
# This script arranges the data into multiple files
#
##


# Libraries

library(dplyr)
library(glue)


# Load the targets

targetsDF <- read.table(
    file = "resources/targets/targets_pw_20.01.22",
    header = T,
    sep = "\t",
    stringsAsFactors = F
)


# Load lm results

chromosomes <- unique(targetsDF$chr)

allResults <- NULL

for (chr in chromosomes) {
    
    print(glue("{Sys.time()}    Importing chromosome {chr}"))
    
    chrResults <- read.table(
        file = glue("C:\\Projects\\placenta_weight\\trio\\moba\\{chr}_20-01-22.lm.gz"),
        header = T,
        sep = "\t",
        stringsAsFactors = F
    )
    
    allResults <- rbind(allResults, chrResults)
    
}

allResults <- allResults %>% distinct()


# Export

targetResults <- allResults %>% 
    filter(
        rsid %in% targetsDF$snp | variantId %in% targetsDF$snp
    )

for (pheno in unique(targetResults$phenotype)) {
    
    print(glue("{Sys.time()}    Exporting {pheno}"))
    
    phenoResults <- targetResults %>% 
        filter(
            phenotype == pheno
        ) %>% 
        select(
            -phenotype
        )
    
    write.table(
        x = phenoResults,
        file = gzfile(glue("docs/results/trio_moba/targets_pw_{pheno}_20.01.22.lm.gz")),
        col.names = T,
        row.names = F,
        quote = F,
        sep = "\t"
    )
    
}

for (i in 1:nrow(targetsDF)) {
    
    rsid <- targetsDF$snp[i]
    chr <- targetsDF$chr[i]
    pos <- targetsDF$start[i]
    
    print(glue("{Sys.time()}    Exporting {rsid}"))
    
    snpResults <- allResults %>% 
        filter(
            contig == chr & abs(position - pos) <= 500000
        )
    
    if (nrow(snpResults) > 0) {
        
        for (pheno in unique(snpResults$phenotype)) {
            
            phenoResults <- snpResults %>% 
                filter(
                    phenotype == pheno
                ) %>% 
                select(
                    -phenotype
                )
            
            write.table(
                x = phenoResults,
                file = gzfile(glue("docs/results/trio_moba/region/targets_pw_{rsid}_{pheno}_20.01.22.lm.gz")),
                col.names = T,
                row.names = F,
                quote = F,
                sep = "\t"
            )
            
        }
    }
}
