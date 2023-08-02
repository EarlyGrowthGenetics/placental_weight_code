#\---------------------------------------------------- Description ----------------------------------------------------#
# Filters and formats the data to fit tensorQTL requirements.
#/---------------------------------------------------------------------------------------------------------------------#

#\----------------------------------------------------- Libraries -----------------------------------------------------#
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(data.table)
#/---------------------------------------------------------------------------------------------------------------------#

#\----------------------------------------------------- Load data -----------------------------------------------------#
args = commandArgs(trailingOnly=TRUE)

load(args[1]) # Load betas
samples <- read.table(args[2])
logFile <- args[3]
savePath <- args[4]
samplesSavePath <- snakemake@output[["samples_out"]]

file.create(logFile)
sink(logFile, append = T, split=T)
#/---------------------------------------------------------------------------------------------------------------------#

#\------------------------------------------------ Create working data ------------------------------------------------#
cat("########################################\nFiltering samples to match WGS data\n########################################\n")
betas <- betas_filtered_var[, colnames(betas_filtered_var) %in% samples$V1]
betas <- betas[,order(as.numeric(colnames(betas)))]
dim(betas)

cat("########################################\nCreating ratio set and annotation\n########################################\n")
RS <- RatioSet(B = betas, annotation=c(array= "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19"))
IlluminaAnnot.EPIC <- getAnnotation(RS)
#/---------------------------------------------------------------------------------------------------------------------#

#\-------------------------------------------------- Temp fix for NAs -------------------------------------------------#
cat("########################################\nDescribe missing data\n########################################\n")
print(paste("Total number of missing values : ", sum(is.na(betas))))
print(paste("Number of CpGs with at least one missing value : ", sum(apply(betas, 1, function(x) any(is.na(x))))))

cat("########################################\nImputing missing data with mean\n########################################\n")
imputeWmean <- function(x){
    x[is.na(x)] <- mean(x, na.rm=T)
    x
}
betas_imputed <- t(apply(betas, 1, imputeWmean))
print(paste("Total number of missing values after imputation : ", sum(is.na(betas_imputed))))
#/---------------------------------------------------------------------------------------------------------------------#

#\-------------------------------------------------- Format into BED --------------------------------------------------#
cat("########################################\nFormat into BED\n########################################\n")
DNAm_bed <- IlluminaAnnot.EPIC[, 1:2]
colnames(DNAm_bed)[2] <- "start"
DNAm_bed$end <- DNAm_bed$start
DNAm_bed$start <- DNAm_bed$start-1
DNAm_bed$ID <- rownames(DNAm_bed)

DNAm_bed <- merge(DNAm_bed, betas_imputed, by=0)
DNAm_bed <- DNAm_bed[, -c(1)]
colnames(DNAm_bed)[1] <- "#chr"
colnames(DNAm_bed)[4] <- "phenotype_id"
#/---------------------------------------------------------------------------------------------------------------------#

#\-------------------------------------------------------- Save -------------------------------------------------------#
fwrite(DNAm_bed, file=savePath, quote=F, sep="\t", row.names=F, nThread=5)
sink()
#/---------------------------------------------------------------------------------------------------------------------#