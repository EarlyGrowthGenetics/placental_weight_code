#\---------------------------------------------------- Description ----------------------------------------------------#
# Filters CpG based on Illumina annotations and a list of cross reactive probes. 
#/---------------------------------------------------------------------------------------------------------------------#

#\----------------------------------------------------- Libraries -----------------------------------------------------#
library(minfi)
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#/---------------------------------------------------------------------------------------------------------------------#

#\----------------------------------------------------- Load data -----------------------------------------------------#
args = commandArgs(trailingOnly=TRUE)

load(args[1]) # Load betas
xReact <- read.csv(args[2])
logFile <- args[3]
savePath <- args[4]

file.create(logFile)
sink(logFile, append = T, split=T)
#/---------------------------------------------------------------------------------------------------------------------#

#\------------------------------------------------ Create working data ------------------------------------------------#
cat("########################################\nCreating ratio set and annotation\n########################################\n")
RS <- RatioSet(B = betas, annotation=c(array= "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19"))
IlluminaAnnot.EPIC <- getAnnotation(RS)
dim(IlluminaAnnot.EPIC)
#/---------------------------------------------------------------------------------------------------------------------#

#\------------------------------------------------------- Filter ------------------------------------------------------#
cat("########################################\nExclude sex chromosomes\n########################################\n")
'%!in%' <- function(x,y)!('%in%'(x,y))
table(IlluminaAnnot.EPIC$chr)
IlluminaAnnot.EPIC <- subset(IlluminaAnnot.EPIC, chr %!in% c("chrX","chrY"))
table(IlluminaAnnot.EPIC$chr)
dim(IlluminaAnnot.EPIC)

cat("########################################\nExclude non-CpG\n########################################\n")
IlluminaAnnot.EPIC <- subset(IlluminaAnnot.EPIC, substring(rownames(IlluminaAnnot.EPIC),1,2) %!in% c("ch", "rs"))
dim(IlluminaAnnot.EPIC)

cat("########################################\nExclude SNP at SBE\n########################################\n")
IlluminaAnnot.EPIC <- subset(IlluminaAnnot.EPIC, is.na(SBE_maf) | SBE_maf < 0.05)
dim(IlluminaAnnot.EPIC)

cat("########################################\nExclude SNP at CpG\n########################################\n")
IlluminaAnnot.EPIC <- subset(IlluminaAnnot.EPIC, is.na(CpG_maf) | CpG_maf < 0.05)
dim(IlluminaAnnot.EPIC)

cat("########################################\nExclude cross-reactive probes\n########################################\n")
nrow(xReact)
rownames(xReact)<-xReact$X
IlluminaAnnot.EPIC <- subset(IlluminaAnnot.EPIC, rownames(IlluminaAnnot.EPIC) %!in%  rownames(xReact))
dim(IlluminaAnnot.EPIC)

cat("########################################\nFilter beta values set\n########################################\n")
betas_filtered <- subset(betas, rownames(betas) %in% rownames(IlluminaAnnot.EPIC))
dim(betas_filtered)
#/---------------------------------------------------------------------------------------------------------------------#

#\-------------------------------------------------------- Save -------------------------------------------------------#
if(!dir.exists(dirname(savePath))) dir.create(dirname(savePath), recursive=T)
save(betas_filtered, file=savePath)
sink()
#/---------------------------------------------------------------------------------------------------------------------#