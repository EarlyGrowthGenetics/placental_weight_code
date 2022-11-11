
library(data.table)
library(dplyr)
library(R.utils)

df <- fread(file = snakemake@input[[1]])

df1 <- fread(file = snakemake@input[[2]])

snplist <- fread(file = snakemake@input[[3]])

snplist <- subset(snplist, select="SNP")

### Maternal Warrington file has rsid already labelled as SNP - need to include Markername =  SNP and SNP = rsid for the fetal file

df1 <- df1 %>% rename (
  Markername = SNPID,
  SNP = RSID,
  A1 = ea,
  A2 = nea,
  FRQ = eaf,
  BETA = beta,
  P = p,
  N = n_offBW
)

##### change A1 and A2 to lower case

df1$A1 <- tolower(df1$A1)
df1$A2 <- tolower(df1$A2)

##### Filter out unnecessary columns

df <- df %>%  select(SNP, A1, A2, FRQ, BETA, P, N)

df1 <- df1 %>% select(SNP, A1, A2, FRQ, BETA, P, N)


### Merge with Hapmap snps to make it faster (alleles that aren't aligned will be corrected in the ldsc code)

df <- merge(snplist, df, by="SNP", all=F)
df1 <- merge(snplist, df1, by="SNP", all=F)

write.table(df, file = snakemake@params[[1]], sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)

gzip(file = snakemake@params[[1]], ext="gz", FUN=gzfile)

write.table(df1, file = snakemake@params[[2]], sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)

gzip(file = snakemake@params[[2]], ext="gz", FUN=gzfile)
