

library(data.table)
library(dplyr)


pwt <- read.table(file = snakemake@input[[1]],  sep="\t", header = TRUE)
ng <- read.table(file = snakemake@input[[2]], sep=",", header = TRUE)


pwt1 <- merge(pwt, ng, by.x="MarkerName", by.y="rsid", all=TRUE)

pwt1 <- pwt1 %>% rename (
    nearestGene = nearest_gene,
    pvalue = P.value,
    CHR = chr,
    POS = pos
)


write.table(pwt1, snakemake@output[[1]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


