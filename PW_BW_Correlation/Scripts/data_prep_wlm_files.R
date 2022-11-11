

library(data.table)
library(dplyr)
library(tidyverse)
library(R.utils)

#### Separate the WLM into fetal, maternal and paternal data frames

wlm <- fread("meta/wlm_sex_gest.gz")

fetga <- fread("meta/pw_fetal_sex_gest.gz")

snplist <- fread('Programs/ldsc/w_hm3.snplist')

snplist <- subset(snplist, select="SNP")


rsid <- subset(fetga, select=c(MarkerName, chr, pos, Allele1, Allele2))

##### Order alleles alphabetically with Allele1 having smallest alphabetical value

rsid <- data.frame(t(apply(rsid[,4:5],1,sort))) %>% 
  rename(Allele1=X1, Allele2=X2) %>% dplyr::bind_cols(rsid[,-4:-5])


rsid <- rsid %>%
  dplyr::rename(rsid = MarkerName)

wlmid <- subset(wlm, select=c(MarkerName))
wlmid$marker <- wlmid$MarkerName

##### Separate MarkerName into chr pos a1 a2

wlmid <- wlmid %>% separate(marker, c("chr", "pos", "Allele1", "Allele2"))

wlmid$Allele1 <- tolower(wlmid$Allele1)
wlmid$Allele2 <- tolower(wlmid$Allele2)

wlmid <- data.frame(t(apply(wlmid[,4:5],1,sort))) %>% 
  rename(Allele1=X1, Allele2=X2) %>% dplyr::bind_cols(wlmid[,-4:-5])

##### Merge wlmid to rsid

wlmid$pos <- as.integer(wlmid$pos)

rsid_wlm <- merge(rsid, wlmid, by=c("chr", "pos", "Allele1", "Allele2"), all.y=TRUE)

rsid_wlm <- rsid_wlm[!duplicated(rsid_wlm$MarkerName),]

rsid_wlm <- subset(rsid_wlm, select=c(MarkerName, rsid))

wlm <- merge(wlm, rsid_wlm, by=c("MarkerName"))

wlmc <- subset(wlm, select=c(rsid, Effect_allele, Other_allele, n_fetal, wlm_beta_fetal, wlm_se_fetal, eaf_fetal))

wlmm <- subset(wlm, select=c(rsid, Effect_allele, Other_allele, n_maternal, wlm_beta_maternal, wlm_se_maternal, eaf_maternal))

wlmd <- subset(wlm, select=c(rsid, Effect_allele, Other_allele, n_paternal, wlm_beta_paternal, wlm_se_paternal, eaf_paternal))



wlmc <- wlmc %>%
  dplyr::rename(
    SNP = rsid,
    A1 = Effect_allele,
    A2 = Other_allele,
    N = n_fetal,
    BETA = wlm_beta_fetal,
    se = wlm_se_fetal,
    FRQ = eaf_fetal
  )


wlmm <- wlmm %>%
  dplyr::rename(
    SNP = rsid,
    A1 = Effect_allele,
    A2 = Other_allele,
    N = n_maternal,
    BETA = wlm_beta_maternal,
    se = wlm_se_maternal,
    FRQ = eaf_maternal
     )


wlmd <- wlmd %>%
  dplyr::rename(
    SNP = rsid,
    A1 = Effect_allele,
    A2 = Other_allele,
    N = n_paternal,
    BETA = wlm_beta_paternal,
    se = wlm_se_paternal,
    FRQ = eaf_paternal
  )

wlmc$Zscore <- wlmc$BETA / wlmc$se

wlmm$Zscore <- wlmm$BETA / wlmm$se

wlmd$Zscore <- wlmd$BETA / wlmd$se

#### calculate the p values for wlm

wlmc$t <- (wlmc$Zscore)^2

wlmc$P <- pchisq(wlmc$t,df=1,lower=F)

wlmm$t <- (wlmm$Zscore)^2

wlmm$P <- pchisq(wlmm$t,df=1,lower=F)

wlmd$t <- (wlmd$Zscore)^2

wlmd$P <- pchisq(wlmd$t,df=1,lower=F)


wlmc <- subset(wlmc, select = c("SNP", "A1", "A2", "BETA", "N", "P", "FRQ"))

wlmm <- subset(wlmm, select = c("SNP", "A1", "A2", "BETA", "N", "P", "FRQ"))

wlmd <- subset(wlmd, select = c("SNP", "A1", "A2", "BETA", "N", "P", "FRQ"))


### Merge with Hapmap snps to make it faster (alleles that aren't aligned will be corrected in the ldsc code)

wlmc <- merge(snplist, wlmc, by="SNP", all=F)

wlmm <- merge(snplist, wlmm, by="SNP", all=F)

wlmd <- merge(snplist, wlmd, by="SNP", all=F)


write.table(wlmc, file = 'Fetal/wlm/fetal_wlm.txt',
            col.names = T, row.names = F, quote = F, sep='\t')

gzip(file = 'Fetal/wlm/fetal_wlm.txt', ext="gz", FUN=gzfile)


write.table(wlmm, file = 'Maternal/wlm/maternal_wlm.txt',
            col.names = T, row.names = F, quote = F, sep='\t')

gzip(file = 'Maternal/wlm/maternal_wlm.txt', ext="gz", FUN=gzfile)


write.table(wlmd, file = 'Paternal/wlm/paternal_wlm.txt',
            col.names = T, row.names = F, quote = F, sep='\t')

gzip(file = 'Paternal/wlm/paternal_wlm.txt', ext="gz", FUN=gzfile)
