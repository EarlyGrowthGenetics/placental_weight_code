#\---------------------------------------------------- Description ----------------------------------------------------#
# Infers population structure on the Gen3G children with the GENESIS pipeline
#/---------------------------------------------------------------------------------------------------------------------#

#\----------------------------------------------------- Libraries -----------------------------------------------------#
library(SNPRelate)
library(GWASTools)
library(GENESIS)
library(ggplot2)
library(ggrepel)
library(scales)
library(stringr)
library(purrr)
library(dplyr)
#/---------------------------------------------------------------------------------------------------------------------#

#\------------------------------------------------ Load data and params -----------------------------------------------#
args = commandArgs(trailingOnly=TRUE)

savePath <- args[1]
inputGDS <- args[2]
outliers <- args[3]
logFile <- args[4]
geneticPCsFile <- args[5]
related_file <- args[6]

file.create(logFile)
sink(logFile, append = T, split=T)

cat("########################################\nLoading data\n########################################\n")
# Create a GdsGenotypeReader which allows to read a GDS file o extract genotype data
geno <- GdsGenotypeReader(filename = inputGDS)
# Create a GenotypeData object from a GdsGenotypeReader connection to a GDS. Contains and allow access to genotype info stored in an R object
# NOTE : this file will have to be closed afterwards
genoData <- GenotypeData(geno)

# Create a SNPGDSFileClass object which provides access to the data stored in a GDS file for methods the directly access the GDS - SNPRelate-specific
# allow.duplicate=T is necessary because a connection to the file has already been opened with GdsGenotypeReader
# NOTE : the file will have to be closed afterwards
gds <- snpgdsOpen(inputGDS, allow.duplicate=T)

# Reading outliers to remove
out_samples <- read.table(outliers)

# Get the list of samples
sample.id <- getScanID(genoData)
thousandGenomID <- sample.id[!is.na(str_extract(sample.id, "^H|N"))]
Gen3GID <- sample.id[!(sample.id %in% thousandGenomID)]
cat(paste0("Removing samples based on ancestry :\n"))
out_samples$V1
Gen3GID <- Gen3GID[!(Gen3GID %in% rownames(out_samples))]
cat(paste0("\nNumber of Gen3G samples left : ", length(Gen3GID), "\n"))
#/---------------------------------------------------------------------------------------------------------------------#

#\------------------------------------------- KING-robust kinship estimation ------------------------------------------#
cat("########################################\nKING-robust kinship estimation\n########################################\n")
dir.create(file.path(savePath,"FirstIter"))

# Prune the SNPs. Using the suggested parameters.
snpset <- snpgdsLDpruning(gds, sample.id=Gen3GID, method="corr", slide.max.bp=10e6, ld.threshold=sqrt(0.1), maf=0.05, missing.rate=0.01, num.thread=10)
cat(paste0("Number of SNPs : ", length(unlist(snpset)), "\n"))

# Infer relatedness with the KING-robust method inplemented in SNPRelate
kingRobust <- snpgdsIBDKING(gds, sample.id=Gen3GID, snp.id=unlist(snpset, use.names=F), num.thread=10)

# Extract the N x N kinship matrix that will be passed to PC-AiR
kingKinship <- kingRobust$kinship
colnames(kingKinship) <- rownames(kingKinship) <- kingRobust$sample.id

# Extract the pairwaise IBD estimates that will allow to explore the data
kingIBD <- snpgdsIBDSelection(kingRobust)

# Plot the kinship
ggplot(kingIBD,aes(x=1:length(ID1),y=kinship))+geom_point()
ggsave(file.path(savePath,"FirstIter","kinship.jpg"))

save(snpset,kingRobust,kingKinship,kingIBD,file=file.path(savePath,"FirstIter","snpset_king.RData"))
#/---------------------------------------------------------------------------------------------------------------------#

#\------------------------------------------------------ PC-AiR -------------------------------------------------------#
cat("########################################\nPC-AiR\n########################################\n")
PcairObjSavePath <- file.path(savePath,"FirstIter","pcair.RData")

# Run PC-AiR, population structure accounting for relatedness
pcairFirstIter <- pcair(gdsobj=genoData, kinobj=kingKinship, divobj=kingKinship, snp.include=unlist(snpset, use.names=F))
# Extract the PCs
pcsFirstIter <- as.data.frame(pcairFirstIter$vectors)
# Extract the % of variance
varProp <- as.data.frame(list(Var=pcairFirstIter$varprop[1:32]))
varProp$Order <- 1:32

# Plot the PCs from the first iteration
ggplot(varProp, aes(x=Order,y=Var*100)) + geom_point() + labs(x="",y="Variance proportion (%)")
ggsave(file.path(savePath,"FirstIter","pcs_var_plot.jpg"))

ggplot(pcsFirstIter,mapping=aes(V1,V2)) + geom_point() + labs(x=paste0("PC1 (",round(varProp[1,1]*100, 3),"%)"),y=paste0("PC2 (",round(varProp[2,1]*100,3),"%)"))
ggsave(file.path(savePath,"FirstIter","pcair_pc1_pc2.jpg"))

ggplot(pcsFirstIter,mapping=aes(V1,V3)) + geom_point() + labs(x=paste0("PC1 (",round(varProp[1,1]*100, 3),"%)"),y=paste0("PC3 (",round(varProp[3,1]*100,3),"%)"))
ggsave(file.path(savePath,"FirstIter","pcair_pc1_pc3.jpg"))

ggplot(pcsFirstIter,mapping=aes(V2,V3)) + geom_point() + labs(x=paste0("PC2 (",round(varProp[2,1]*100, 3),"%)"),y=paste0("PC3 (",round(varProp[3,1]*100,3),"%)"))
ggsave(file.path(savePath,"FirstIter","pcair_pc2_pc3.jpg"))

save(pcairFirstIter,pcsFirstIter,varProp,file=PcairObjSavePath)
#/---------------------------------------------------------------------------------------------------------------------#

#\----------------------------------------------------- PC-Relate -----------------------------------------------------#
cat("########################################\nPC-Relate\n########################################\n")
PcrelateObjSavePath <- file.path(savePath,"SecondIter","pcrelate.RData")
dir.create(file.path(savePath,"SecondIter"))

# Create a GenomicDataIterator reading only the pruned SNPs
genoDataIter <- GenotypeBlockIterator(genoData, snpInclude=unlist(snpset, use.names=F))
# PC-Realte, Compute relatedness accounting for ancestry (first 4 PC of PC-AiR)
pcrel <- pcrelate(genoDataIter, sample.include=Gen3GID, pcs=pcairFirstIter$vectors[,1:4], training.set=pcairFirstIter$unrels)
# Extract a GRM, N x N kinship matrix to use in the second iteration of PC-AiR
pcrelGRM <- pcrelateToMatrix(pcrel)
# Plot the relatedness
pcrel$kinBtwn$ID1 <- sub("Cordblood", "ONT", pcrel$kinBtwn$ID1)
pcrel$kinBtwn$ID2 <- sub("Cordblood", "ONT", pcrel$kinBtwn$ID2)
ggplot(pcrel$kinBtwn,aes(x=1:length(ID1),y=kin,label=ifelse(kin>0.1,paste(ID1,ID2,sep=" - "),"")))+geom_point()+geom_text_repel(force=3,box.padding=2, max.overlaps = Inf)+labs(x="",y="kinship")+ylim(c(-0.05,0.6))
ggsave(file.path(savePath,"SecondIter","kinship.jpg"))
ggplot(pcrel$kinBtwn,aes(x=k0,y=kin,label=ifelse(kin>0.1,paste(ID1,ID2,sep=" - "),"")))+geom_point()+geom_text_repel(force=7, max.overlaps = Inf, xlim=c(0.6, Inf))
ggsave(file.path(savePath,"SecondIter","kinship_IBD0.jpg"))
ggplot(pcrel$kinBtwn,aes(x=k0,y=k2,label=ifelse(k0<0.6,paste(ID1,ID2,sep=" - "),"")))+geom_point()+geom_text_repel(force=7, max.overlaps = Inf, xlim=c(0.6, Inf))
ggsave(file.path(savePath,"SecondIter","IBD0_IBD2.jpg"))

# Save relatedness table
write.table(pcrel$kinBtwn, file=file.path(savePath,"SecondIter","kinship.tsv"))
save(pcrel,pcrelGRM,file=PcrelateObjSavePath)
#/---------------------------------------------------------------------------------------------------------------------#

#\------------------------------------------------- PC-AiR - Second iteration -----------------------------------------#
cat("########################################\nPC-AiR - Second iteration\n########################################\n")
PcairSecItObjSavePath <- file.path(savePath,"SecondIter","pcair.RData")

# PC-AiR
pcairSecondIter <- pcair(gdsobj=genoData, kinobj=pcrelGRM, divobj=pcrelGRM, snp.include=unlist(snpset, use.names=F))
pcsSecondIter <- as.data.frame(pcairSecondIter$vectors)
varPropSecondIter <- as.data.frame(list(Var=pcairSecondIter$varprop[1:32]))
varPropSecondIter$Order <- 1:32

# Plot the PCs
ggplot(varPropSecondIter, aes(x=Order,y=Var*100)) + geom_point() + labs(x="",y="Variance proportion (%)") + scale_y_continuous(breaks = round(seq(min(varPropSecondIter$Var*100), max(varPropSecondIter$Var*100), length.out = 5), 3))
ggsave(file.path(savePath,"SecondIter","pcs_var_plot.jpg"))

pcsSecondIter$ID <- rownames(pcsSecondIter)

# Save the first PCs
write.table(pcsSecondIter[, c("ID", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10")], file=geneticPCsFile, row.names=F)

# Calculate threshold for outliers
out_PC1 <- c(mean(pcsSecondIter$V1)-3*sd(pcsSecondIter$V1), mean(pcsSecondIter$V1)+3*sd(pcsSecondIter$V1))
out_PC2 <- c(mean(pcsSecondIter$V2)-3*sd(pcsSecondIter$V2), mean(pcsSecondIter$V2)+3*sd(pcsSecondIter$V2))
out_PC3 <- c(mean(pcsSecondIter$V3)-3*sd(pcsSecondIter$V3), mean(pcsSecondIter$V3)+3*sd(pcsSecondIter$V3))

ggplot(pcsSecondIter,mapping=aes(V1,V2,label=ifelse(between(V1, out_PC1[1], out_PC1[2]) & between(V2, out_PC2[1], out_PC2[2]),"",ID))) + geom_point() + labs(x=paste0("PC1 (",round(varPropSecondIter[1,1]*100, 3),"%)"),y=paste0("PC2 (",round(varPropSecondIter[2,1]*100,3),"%)")) + geom_text_repel()
ggsave(file.path(savePath,"SecondIter","pcair_pc1_pc2.jpg"))

ggplot(pcsSecondIter,mapping=aes(V1,V3,label=ifelse(between(V1, out_PC1[1], out_PC1[2]) & between(V3, out_PC3[1], out_PC3[2]),"",ID))) + geom_point() + labs(x=paste0("PC1 (",round(varPropSecondIter[1,1]*100, 3),"%)"),y=paste0("PC3 (",round(varPropSecondIter[3,1]*100,3),"%)")) + geom_text_repel()
ggsave(file.path(savePath,"SecondIter","pcair_pc1_pc3.jpg"))

ggplot(pcsSecondIter,mapping=aes(V2,V3,label=ifelse(between(V2, out_PC2[1], out_PC2[2]) & between(V3, out_PC3[1], out_PC3[2]),"",ID))) + geom_point() + labs(x=paste0("PC2 (",round(varPropSecondIter[2,1]*100, 3),"%)"),y=paste0("PC3 (",round(varPropSecondIter[3,1]*100,3),"%)")) + geom_text_repel()
ggsave(file.path(savePath,"SecondIter","pcair_pc2_pc3.jpg"))

save(pcairSecondIter,pcsSecondIter,varPropSecondIter,file=PcairSecItObjSavePath)
#/---------------------------------------------------------------------------------------------------------------------#

#\---------------------------------------------------- Find related ---------------------------------------------------#
cat("########################################\nFinding related pairs\n########################################\n")
dir.create(file.path(savePath,"Related"))
# Create scale for relatedness
# 2^-5.5 ~=0.0221, 2^-3.5 ~=0.0884, 2^-2.5~=0.177, 2^-1.5~=0.354
rel_scale <- scale_color_manual(name = "Estimated rel", values = c("(-Inf,0.0221]"=hue_pal()(5)[1], "(0.0221,0.0884]"=hue_pal()(5)[2], "(0.0884,0.177]"=hue_pal()(5)[3], "(0.177,0.354]"=hue_pal()(5)[4], "(0.354, Inf]"=hue_pal()(5)[5]), labels = c("< second cousins", "first cousins", "half-sibling", "siblings", "self"))

ggplot(pcrel$kinBtwn,aes(x=1:length(ID1),y=kin))+geom_point(aes(colour=cut(kin, c(-Inf, 2^-5.5, 2^-3.5, 2^-2.5, 2^-1.5, Inf))))+labs(x="",y="kinship")+geom_hline(yintercept=c(2^-5.5)) + rel_scale
ggsave(file.path(savePath,"Related","kinship_coloured.jpg"))

related <- subset(pcrel$kinBtwn, kin > 2^-5.5)
write.table(related, file=related_file, row.names=F)
#/---------------------------------------------------------------------------------------------------------------------#

#\---------------------------------------------------- Close files ----------------------------------------------------#
snpgdsClose(gds)
close(genoData)
sink()
#/---------------------------------------------------------------------------------------------------------------------#