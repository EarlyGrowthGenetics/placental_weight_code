
library(data.table)
library(dplyr)
library(ggplot2)

# first use head and grep to extract top hits from cleaned individual cohorts
# zcat EGG_HRC_BW6.PW.child.sex_gest.Roskilde.european.SES.20200403.txt.gz|head -n 1 > Roskilde.txt
# zgrep -f top.csv EGG_HRC_BW6.PW.child.sex_gest.Roskilde.european.SES.20200403.txt.gz >> Roskilde.txt

# load paternal meta summary stats

Meta <- fread('/child_gest/LDSC_Files/Paternal/Clean/pw_paternal_sex_gest.gz')

# create marker column

Meta$marker <- paste(Meta$chr, Meta$pos, sep=':')

# rename columns

Meta <- Meta %>% rename(
                    SNP = MarkerName,
                    A1 = Allele1,
                    A2 = Allele2,
                    FRQ = Freq1,
                    BETA = Effect,
                    SE = StdErr,
                    P = 'P-value',
                    N = TOTALSAMPLESIZE
    ) 

# load top hits and filter out chr 22 where needed (has the same position and identified using grep merge will double check every variant)

th <- fread('/child_gest/LDSC_Files/paternal_top.csv', header=F)

th <- th %>% rename(marker = V1)

# keep only two top variants from paternal Meta

Meta <- inner_join(Meta, th, by = 'marker')

MEta <- Meta %>% select("SNP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N")

# Set the path to the folder containing the .txt files
path <- '/child_gest/LDSC_Files/Paternal/top_hits_files'

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
MoBa$marker <- paste(MoBa$CHR, MoBa$POS, sep=':')
DBDS$marker <- paste(DBDS$CHR, DBDS$POS, sep=':')
EFSOCH <- EFSOCH %>%  rename(marker = SNPID)
GDAffy <- GDAffy %>%  rename(marker = SNPID)
GDIllum <- GDIllum %>%  rename(marker = SNPID)
GSAV2 <- GSAV2 %>%  rename(marker = SNPID)
HUNT$marker <- paste(HUNT$CHR, HUNT$POS, sep=':')
RaineStudy$marker <- paste(RaineStudy$CHR, RaineStudy$POS, sep=':')


# create a list of dataframes

df_list <- list(DBDS, EFSOCH, GDAffy, GDIllum, GSAV2, HUNT, Meta, MoBa, RaineStudy)


# merge dataframes with th dataframe to check only valid variants included

joined_data <- lapply(df_list, function(x) inner_join(th, x, by='marker', all=F))

names(joined_data) <- c("DBDS", "EFSOCH", "GDAffy", "GDIllum", "GSAV2", "HUNT", "Meta", "MoBa", "RaineStudy")

list2env(joined_data, envir = .GlobalEnv)


cohorts <- list(DBDS, EFSOCH, GDAffy, GDIllum, GSAV2, HUNT, RaineStudy)


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

names(df_list) <- c("DBDS", "EFSOCH", "GDAffy", "GDIllum", "GSAV2", "HUNT", "RaineStudy")

list2env(df_list, envir = .GlobalEnv)

cohorts <- list(DBDS, EFSOCH, GDAffy, GDIllum, GSAV2, HUNT, Meta, MoBa, RaineStudy)

names(cohorts) <- list("DBDS", "EFSOCH", "GDAffy", "GDIllum", "GSAV2", "HUNT", "Meta", "MoBa", "RaineStudy")

var_list  <- c("SNP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N")

# keep only variables needed
df_list <- lapply(cohorts, function(x) x[, var_list, with = FALSE])

names(df_list) <- c("DBDS", "EFSOCH", "GDAffy", "GDIllum", "GSAV2", "HUNT", "Meta", "MoBa", "RaineStudy")

list2env(df_list, envir = .GlobalEnv)

# harmonise all dataframes against the Meta dataframe

data_list <- list(DBDS, EFSOCH, GDAffy, GDIllum, GSAV2, HUNT, Meta, MoBa, RaineStudy)

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

names(data_list_modified) <- c("DBDS", "EFSOCH", "GDAffy", "GDIllum", "GSAV2", "HUNT", "MoBa", "RaineStudy")

list2env(data_list_modified, envir = .GlobalEnv)


df_list <- list(DBDS, EFSOCH, GDAffy, GDIllum, GSAV2, HUNT, Meta, MoBa, RaineStudy)

names(df_list) <- c("DBDS", "EFSOCH", "GDAffy", "GDIllum", "GSAV2", "HUNT", "Meta", "MoBa", "RaineStudy")

# Use lapply to write each dataframe to a .txt file in the specified folder

lapply(seq_along(df_list), function(i){
  # Extract the name of the dataframe
  df_name <- names(df_list)[i]
  # Write the dataframe to a .csv file in the specified folder
  write.table(df_list[[i]], file = paste0('/child_gest/LDSC_Files/Paternal/Cleaned_data', '/', df_name, '.txt'),sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
})


