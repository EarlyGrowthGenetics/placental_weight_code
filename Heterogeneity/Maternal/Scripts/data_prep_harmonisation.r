
library(data.table)
library(dplyr)
library(ggplot2)

# first use head and grep to extract top hits from cleaned individual cohorts
# zcat EGG_HRC_BW6.PW.child.sex_gest.Roskilde.european.SES.20200403.txt.gz|head -n 1 > Roskilde.txt
# zgrep -f top.csv EGG_HRC_BW6.PW.child.sex_gest.Roskilde.european.SES.20200403.txt.gz >> Roskilde.txt

# load top hits and filter out chr 22 where needed (has the same position and identified using grep merge will double check every variant)

th <- fread('/child_gest/LDSC_Files/maternal_top.csv', header=F)

# Set the path to the folder containing the .txt files
path <- '/child_gest/LDSC_Files/Maternal/top_hits_files'

# Get a list of all .txt files in the folder
txt_files <- list.files(path, pattern = "\\.txt$", full.names = TRUE)

# Loop through each .txt file, read it into a dataframe, and assign the file name as the dataframe name
for (file in txt_files) {
  df <- read.table(file, header = TRUE, sep = "\t")
  file_name <- gsub("\\.txt$", "", basename(file))
  assign(file_name, df)
}

rm(df)

# create a marker variable for each cohort
th <- th %>% rename(marker = V1)
MoBa$marker <- paste(MoBa$CHR, MoBa$POS, sep=':')
ALSPAC <- ALSPAC %>% rename(marker = SNPID)
CHB$marker <- paste(CHB$CHR, CHB$POS, sep=':')
DNBC$marker <- paste(DNBC$CHR, DNBC$POS, sep=':')
DBDS$marker <- paste(DBDS$CHR, DBDS$POS, sep=':')
EFSOCH <- EFSOCH %>%  rename(marker = SNPID)
GDAffy <- GDAffy %>%  rename(marker = SNPID)
Gen3G$marker <- paste(Gen3G$CHR, Gen3G$POS, sep=':')
GOYA_obese_children$marker <- paste(GOYA_obese_children$CHR, GOYA_obese_children$POS, sep=':')
GOYA_obese_mothers$marker <- paste(GOYA_obese_mothers$CHR, GOYA_obese_mothers$POS, sep=':')
GOYA_control_mothers$marker <- paste(GOYA_control_mothers$CHR, GOYA_control_mothers$POS, sep=':')
GSAV2 <- GSAV2 %>%  rename(marker = SNPID)
HUNT$marker <- paste(HUNT$CHR, HUNT$POS, sep=':')
Mono_ctrls$marker <- paste(Mono_ctrls$CHR, Mono_ctrls$POS, sep=':')
PE_ctrls$marker <- paste(PE_ctrls$CHR, PE_ctrls$POS, sep=':')
PPD_moms$marker <- paste(PPD_moms$CHR, PPD_moms$POS, sep=':')
RaineStudy$marker <- paste(RaineStudy$CHR, RaineStudy$POS, sep=':')


#GSAV2 has weird decimal placed numbers all should be 195

GSAV2$N <- 168


# create a list of dataframes

df_list <- list(ALSPAC, CHB, DNBC, DBDS, EFSOCH, GDAffy, Gen3G, GOYA_control_mothers, GOYA_obese_children, GOYA_obese_mothers,
                GSAV2, HUNT, Mono_ctrls, PE_ctrls, PPD_moms, Meta, MoBa, RaineStudy)


# merge dataframes with th dataframe to check only valid variants included

joined_data <- lapply(df_list, function(x) inner_join(th, x, by='marker', all=F))

names(joined_data) <- c("ALSPAC", "CHB", "DNBC", "DBDS", "EFSOCH", "GDAffy", "Gen3G", "GOYA_control_mothers", "GOYA_obese_children", "GOYA_obese_mothers", 
                    "GSAV2", "HUNT", "Mono_ctrls", "PE_ctrls", "PPD_moms", "Meta", "MoBa", "RaineStudy")

list2env(joined_data, envir = .GlobalEnv)


cohorts <- list(ALSPAC, CHB, DNBC, DBDS, EFSOCH, GDAffy, Gen3G, GOYA_control_mothers, GOYA_obese_children, GOYA_obese_mothers,
                GSAV2, HUNT, Mono_ctrls, PE_ctrls, PPD_moms, RaineStudy)


var_list <- c("marker", "EFFECT_ALLELE", "NON_EFFECT_ALLELE", "EAF", "BETA", "SE", "PVAL", "N")

MoBa <- MoBa %>% select('marker', 'SNPID', 'EFFECT_ALLELE', 'NON_EFFECT_ALLELE', 'EAF', 'BETA', 'SE', 'PVAL', 'N') %>% 
   rename (
          Markername = marker,
          SNP = SNPID,
          A1 = EFFECT_ALLELE,
          A2 = NON_EFFECT_ALLELE,
          FRQ = EAF,
          P = PVAL
) 


#get rsid from MoBa to merge with other cohorts

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

# change alleles to lowercase in MoBa

MoBa$A1 <- tolower(MoBa$A1)
MoBa$A2 <- tolower(MoBa$A2)

#create names for and change from df_list to dataframes in the global environment

names(df_list) <- c("ALSPAC", "CHB", "DNBC", "DBDS", "EFSOCH", "GDAffy", "Gen3G", "GOYA_control_mothers", "GOYA_obese_children", "GOYA_obese_mothers", 
                    "GSAV2", "HUNT", "Mono_ctrls", "PE_ctrls", "PPD_moms", "RaineStudy")

list2env(df_list, envir = .GlobalEnv)

cohorts <- list(ALSPAC, CHB, DNBC, DBDS, EFSOCH, GDAffy, Gen3G, GOYA_control_mothers, GOYA_obese_children, GOYA_obese_mothers,
                GSAV2, HUNT, Meta, MoBa, Mono_ctrls, PE_ctrls, PPD_moms, RaineStudy)

names(cohorts) <- list("ALSPAC", "CHB", "DNBC", "DBDS", "EFSOCH", "GDAffy", "Gen3G", "GOYA_control_mothers", "GOYA_obese_children", "GOYA_obese_mothers", 
                       "GSAV2", "HUNT", "Meta", "MoBa", "Mono_ctrls", "PE_ctrls", "PPD_moms", "RaineStudy")

var_list  <- c("SNP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N")

# keep only variables needed
df_list <- lapply(cohorts, function(x) x[, var_list, with = FALSE])

names(df_list) <- c("ALSPAC", "CHB", "DNBC", "DBDS", "EFSOCH", "GDAffy", "Gen3G", "GOYA_control_mothers", "GOYA_obese_children", "GOYA_obese_mothers", 
                    "GSAV2", "HUNT", "Meta", "MoBa", "Mono_ctrls", "PE_ctrls", "PPD_moms", "RaineStudy")

list2env(df_list, envir = .GlobalEnv)

# harmonise all dataframes against the Meta dataframe

data_list <- list(ALSPAC, CHB, DNBC, DBDS, EFSOCH, GDAffy, Gen3G, GOYA_control_mothers, GOYA_obese_children, GOYA_obese_mothers,
                  GSAV2, HUNT, MoBa, Mono_ctrls, PE_ctrls, PPD_moms, RaineStudy)

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

names(data_list_modified) <- c("ALSPAC", "CHB", "DNBC", "DBDS", "EFSOCH", "GDAffy", "Gen3G", "GOYA_control_mothers", "GOYA_obese_children", "GOYA_obese_mothers", 
                               "GSAV2", "HUNT", "MoBa", "Mono_ctrls", "PE_ctrls", "PPD_moms", "RaineStudy")

list2env(data_list_modified, envir = .GlobalEnv)


df_list <- list(ALSPAC, CHB, DNBC, DBDS, EFSOCH, GDAffy, Gen3G, GOYA_control_mothers, GOYA_obese_children, GOYA_obese_mothers,
                GSAV2, HUNT, Meta, MoBa, Mono_ctrls, PE_ctrls, PPD_moms, RaineStudy)

names(df_list) <- c("ALSPAC", "CHB", "DNBC", "DBDS", "EFSOCH", "GDAffy", "Gen3G", "GOYA_control_mothers", "GOYA_obese_children", "GOYA_obese_mothers", 
                    "GSAV2", "HUNT", "Meta", "MoBa", "Mono_ctrls", "PE_ctrls", "PPD_moms", "RaineStudy")

# Use lapply to write each dataframe to a .txt file in the specified folder

lapply(seq_along(df_list), function(i){
  # Extract the name of the dataframe
  df_name <- names(df_list)[i]
  # Write the dataframe to a .csv file in the specified folder
  write.table(df_list[[i]], file = paste0('/child_gest/LDSC_Files/Maternal/Cleaned_data', '/', df_name, '.txt'),sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
})


