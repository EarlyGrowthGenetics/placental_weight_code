
library(data.table)
library(dplyr)
library(ggplot2)

# first use head and grep to extract top hits from cleaned individual cohorts
# zcat EGG_HRC_BW6.PW.child.sex_gest.Roskilde.european.SES.20200403.txt.gz|head -n 1 > Roskilde.txt
# zgrep -f top.csv EGG_HRC_BW6.PW.child.sex_gest.Roskilde.european.SES.20200403.txt.gz >> Roskilde.txt

# load top hits and filter out chr 22 where needed (has the same position and identified using grep merge will double check every variant)

th <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top.csv', header=F)

Meta <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/Meta.txt')
MoBa <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/MoBa.txt')
ALSPAC <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/ALSPAC.txt')
DNBC <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/DNBC.txt')
EFSOCH <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/EFSOCH.txt')
GENR <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/GENR.txt')
GOYA_obese_children <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/GOYA_obese_children.txt')
GOYA_obese_mothers <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/GOYA_obese_mothers.txt')
GOYA_control_mothers <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/GOYA_control_mothers.txt')
FS <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/FS.txt')
GSAV2 <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/GSAV2.txt')
IHPS_casesCtrls_MEGAex <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/IHPS_casesCtrls_MEGAex.txt')
IHPS_CIDR <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/IHPS_CIDR.txt')
INMAGSA <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/INMAGSA.txt')
INMAOmni <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/INMAOmni.txt')
iPSYCH <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/iPSYCH.txt')
NFBC1966 <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/NFBC1966.txt')
NFBC1986 <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/NFBC1986.txt')
OPI <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/OPI.txt')
PANIC <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/PANIC.txt')
RaineStudy <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/RaineStudy.txt')
Roskilde <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/top_hits_files/Roskilde.txt')


th <- th %>% rename(marker = V1)
Meta <- Meta %>% rename(marker ="chr:pos")
MoBa$marker <- paste(MoBa$CHR, MoBa$POS, sep=':')
ALSPAC <- ALSPAC %>%  rename(marker = SNPID)
DNBC$marker <- paste(DNBC$CHR, DNBC$POS, sep=':')
EFSOCH <- EFSOCH %>%  rename(marker = SNPID)
GENR <- GENR %>%  rename(marker = SNPID)
GOYA_obese_children$marker <- paste(GOYA_obese_children$CHR, GOYA_obese_children$POS, sep=':')
GOYA_obese_mothers$marker <- paste(GOYA_obese_mothers$CHR, GOYA_obese_mothers$POS, sep=':')
GOYA_control_mothers$marker <- paste(GOYA_control_mothers$CHR, GOYA_control_mothers$POS, sep=':')
FS$marker <- paste(FS$CHR, FS$POS, sep=':')
GSAV2 <- GSAV2 %>%  rename(marker = SNPID)
IHPS_casesCtrls_MEGAex$marker <- paste(IHPS_casesCtrls_MEGAex$CHR, IHPS_casesCtrls_MEGAex$POS, sep=':')
IHPS_CIDR$marker <- paste(IHPS_CIDR$CHR, IHPS_CIDR$POS, sep=':')
INMAGSA$marker <- paste(INMAGSA$CHR, INMAGSA$POS, sep=':')
INMAOmni$marker <- paste(INMAOmni$CHR, INMAOmni$POS, sep=':')
iPSYCH$marker <- paste(iPSYCH$CHR, iPSYCH$POS, sep=':')
NFBC1966 <- NFBC1966 %>%  rename(marker = SNPID)
NFBC1986 <- NFBC1986 %>%  rename(marker = SNPID)
OPI$marker <- paste(OPI$CHR, OPI$POS, sep=':')
PANIC$marker <- paste(PANIC$CHR, PANIC$POS, sep=':')
RaineStudy$marker <- paste(RaineStudy$CHR, RaineStudy$POS, sep=':')
Roskilde$marker <- paste(Roskilde$CHR, Roskilde$POS, sep=':')

#GSAV2 has weird decimal placed numbers all should be 195

GSAV2$N <- 195

names(iPSYCH)[2] <- "SNP_Name"

# merge dataframes with th dataframe to check only valid variants included (One in chr 22 had the same position and got picked up with grep)
Meta <- inner_join(th, Meta, by='marker', all=F)
MoBa <- inner_join(th, MoBa, by='marker', all=F)
ALSPAC <- inner_join(th, ALSPAC, by='marker', all=F)
DNBC <- inner_join(th, DNBC, by='marker', all=F)
EFSOCH <- inner_join(th, EFSOCH, by='marker', all=F)
GENR <- inner_join(th, GENR, by='marker', all=F)
GOYA_obese_children <- inner_join(th, GOYA_obese_children, by='marker', all=F)
GOYA_obese_mothers <- inner_join(th, GOYA_obese_mothers, by='marker', all=F)
GOYA_control_mothers <- inner_join(th, GOYA_control_mothers, by='marker', all=F)
FS <- inner_join(th, FS, by='marker', all=F)
GSAV2 <- inner_join(th, GSAV2, by='marker', all=F)
IHPS_casesCtrls_MEGAex <- inner_join(th, IHPS_casesCtrls_MEGAex, by='marker', all=F)
IHPS_CIDR <- inner_join(th, IHPS_CIDR, by='marker', all=F)
INMAGSA <- inner_join(th, INMAGSA, by='marker', all=F)
INMAOmni <- inner_join(th, INMAOmni, by='marker', all=F)
iPSYCH <- inner_join(th, iPSYCH, by='marker', all=F)
NFBC1966 <- inner_join(th, NFBC1966, by='marker', all=F)
NFBC1986 <- inner_join(th, NFBC1986, by='marker', all=F)
OPI <- inner_join(th, OPI, by='marker', all=F)
PANIC <- inner_join(th, PANIC, by='marker', all=F)
RaineStudy <- inner_join(th, RaineStudy, by='marker', all=F)
Roskilde <- inner_join(th, Roskilde, by='marker', all=F)

Meta <- Meta %>%rename (
  Markername = marker,
  SNP = MarkerName,
  A1 = Allele1,
  A2 = Allele2,
  FRQ = Freq1,
  BETA = Effect,
  SE = StdErr,
  P = 'P-value',
  N = TOTALSAMPLESIZE
) %>%  select('Markername', 'SNP', 'A1', 'A2', 'FRQ', 'BETA', 'SE','P', 'N')

cohorts <- list(ALSPAC, DNBC, EFSOCH, FS, GENR, GOYA_control_mothers, GOYA_obese_children,  GOYA_obese_mothers,
                GSAV2, IHPS_casesCtrls_MEGAex, IHPS_CIDR, INMAGSA, INMAOmni, iPSYCH, 
                NFBC1966, NFBC1986, OPI, PANIC, RaineStudy, Roskilde)


var_list <- c("marker", "EFFECT_ALLELE", "NON_EFFECT_ALLELE", "EAF", "BETA", "SE","PVAL", "N")

MoBa <- MoBa %>% select('marker', 'SNPID', 'EFFECT_ALLELE', 'NON_EFFECT_ALLELE', 'EAF', 'BETA', 'SE', 'PVAL', 'N') %>% 
   rename (
          Markername = marker,
          SNP = SNPID,
          A1 = EFFECT_ALLELE,
          A2 = NON_EFFECT_ALLELE,
          FRQ = EAF,
          P = PVAL
) 

#get rsid from MoBa to merge with pther cohorts

snp <- MoBa %>% select(Markername, SNP)

# keep only variables needed
df_list <- lapply(cohorts, function(x) x[, var_list, with = FALSE])

# create a list of new variable names
new_names <- list(c("Markername", "A1", "A2", "FRQ"))

# use lapply to rename variables in each dataframe
df_list <- lapply(df_list, function(x) {
  rename(x,           
         Markername = marker,
         A1 = EFFECT_ALLELE,
         A2 = NON_EFFECT_ALLELE,
         FRQ = EAF,
         P = PVAL)
})

# merge rsids into each cohort. SNP taken from MoBa
df_list <- lapply(df_list, function(x) merge(x, snp, by = "Markername"))


# change alleles to lowercase
df_list <- lapply(df_list, function(x) {
  x$A1 <- tolower(x$A1)
  x$A2 <- tolower(x$A2)
  return(x)
})

MoBa$A1 <- tolower(MoBa$A1)
MoBa$A2 <- tolower(MoBa$A2)

names(df_list) <- c("ALSPAC", "DNBC", "EFSOCH", "FS", "GENR", "GOYA_control_mothers", "GOYA_obese_children", "GOYA_obese_mothers", 
                      "GSAV2", "IHPS_casesCtrls_MEGAex", "IHPS_CIDR", "INMAGSA", "INMAOmni", "iPSYCH", "NFBC1966", "NFBC1986", 
                       "OPI", "PANIC", "RaineStudy", "Roskilde")

list2env(df_list, envir = .GlobalEnv)

cohorts <- list(ALSPAC, DNBC, EFSOCH, FS, GENR, GOYA_control_mothers, GOYA_obese_children,  GOYA_obese_mothers,
                  GSAV2, IHPS_casesCtrls_MEGAex, IHPS_CIDR, INMAGSA, INMAOmni, iPSYCH, Meta, MoBa,
                    NFBC1966, NFBC1986, OPI, PANIC, RaineStudy, Roskilde)

names(cohorts) <- list("ALSPAC", "DNBC", "EFSOCH", "FS", "GENR", "GOYA_control_mothers", "GOYA_obese_children", "GOYA_obese_mothers",
                         "GSAV2", "IHPS_casesCtrls_MEGAex", "IHPS_CIDR", "INMAGSA", "INMAOmni", "iPSYCH", "Meta", "MoBa",
                           "NFBC1966", "NFBC1986", "OPI", "PANIC", "RaineStudy", "Roskilde")

var_list  <- c("SNP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N")

# keep only variables needed
df_list <- lapply(cohorts, function(x) x[, var_list, with = FALSE])

names(df_list) <- c("ALSPAC", "DNBC", "EFSOCH", "FS", "GENR", "GOYA_control_mothers", "GOYA_obese_children", "GOYA_obese_mothers", 
                      "GSAV2", "IHPS_casesCtrls_MEGAex", "IHPS_CIDR", "INMAGSA", "INMAOmni", "iPSYCH", "Meta", "MoBa", "NFBC1966", "NFBC1986", 
                        "OPI", "PANIC", "RaineStudy", "Roskilde")

list2env(df_list, envir = .GlobalEnv)

# harmonise all dataframes against the Meta dataframe

data_list <- list(ALSPAC, DNBC, EFSOCH, FS, GENR, GOYA_control_mothers, GOYA_obese_children, GOYA_obese_mothers,
                    GSAV2, IHPS_casesCtrls_MEGAex, IHPS_CIDR, INMAGSA, INMAOmni, iPSYCH, MoBa,
                      NFBC1966, NFBC1986, OPI, PANIC, RaineStudy, Roskilde)

operation_func <- function(df) {
  # Merge the dataframe with Meta
  merged_df <- merge(Meta, df, by = "SNP", all.y = TRUE)
  
  # Find the rows where A1 in Meta and df do not match
  mismatch_rows <- which(merged_df$A1.x != merged_df$A1.y)
  
  # Swap the values in A1 and A2 for those rows in df
  merged_df$A1.y[mismatch_rows] <- merged_df$A1.x[mismatch_rows]
  merged_df$A2.y[mismatch_rows] <- merged_df$A2.x[mismatch_rows]
  
  # Multiply BETA by -1 for those rows in df
  merged_df$BETA.y[mismatch_rows] <- -1*merged_df$BETA.y[mismatch_rows]
  
  # Need to change the frequency for those mismatched too
  merged_df$FRQ.y[mismatch_rows] <- 1-merged_df$FRQ.y[mismatch_rows]
  
  # Subset the merged dataframe to get the rows that belong to df
  df_subset <- merged_df %>%  select(c("SNP", "A1.y", "A2.y", "FRQ.y", "BETA.y", "SE.y", "P.y", "N.y"))
  
  # Rename the columns of df_subset
  names(df_subset) <- c("SNP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N")
  return(df_subset)
}

# Apply the function to each dataframe in the list
data_list_modified <- lapply(data_list , operation_func)

names(data_list_modified) <- c("ALSPAC", "DNBC", "EFSOCH", "FS", "GENR", "GOYA_control_mothers", "GOYA_obese_children", "GOYA_obese_mothers", 
                                 "GSAV2", "IHPS_casesCtrls_MEGAex", "IHPS_CIDR", "INMAGSA", "INMAOmni", "iPSYCH", "MoBa", "NFBC1966", "NFBC1986", 
                                   "OPI", "PANIC", "RaineStudy", "Roskilde")

list2env(data_list_modified, envir = .GlobalEnv)


df_list <- list(ALSPAC, DNBC, EFSOCH, FS, GENR, GOYA_control_mothers, GOYA_obese_children, GOYA_obese_mothers,
                  GSAV2, IHPS_casesCtrls_MEGAex, IHPS_CIDR, INMAGSA, INMAOmni, iPSYCH, Meta, MoBa,
                    NFBC1966, NFBC1986, OPI, PANIC, RaineStudy, Roskilde)

names(df_list) <- c("ALSPAC", "DNBC", "EFSOCH", "FS", "GENR", "GOYA_control_mothers", "GOYA_obese_children", "GOYA_obese_mothers", 
                      "GSAV2", "IHPS_casesCtrls_MEGAex", "IHPS_CIDR", "INMAGSA", "INMAOmni", "iPSYCH", "Meta", "MoBa", "NFBC1966", "NFBC1986", 
                        "OPI", "PANIC", "RaineStudy", "Roskilde")

# Use lapply to write each dataframe to a .txt file in the specified folder

lapply(seq_along(df_list), function(i){
  # Extract the name of the dataframe
  df_name <- names(df_list)[i]
  # Write the dataframe to a .csv file in the specified folder
  write.table(df_list[[i]], file = paste0('/home/christopher/Desktop/child_gest/LDSC_Files/cleaned_data', '/', df_name, '.txt'),sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
})


