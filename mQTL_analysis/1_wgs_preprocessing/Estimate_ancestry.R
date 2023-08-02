#\---------------------------------------------------- Description ----------------------------------------------------#
# Infers population structure on the Gen3G children with the GENESIS pipeline based on the 1000G dataset
#/---------------------------------------------------------------------------------------------------------------------#

#\----------------------------------------------------- Libraries -----------------------------------------------------#
library(SNPRelate)
library(GWASTools)
library(GENESIS)
library(ggplot2)
library(stringr)
library(purrr)
#/---------------------------------------------------------------------------------------------------------------------#

#\------------------------------------------------ Load data and params -----------------------------------------------#
args = commandArgs(trailingOnly=TRUE)
savePath <- args[1]
inputGDS <- args[2]
logFile <- args[3]
metaDataFile <- args[4]
outliers_file <- args[5]

file.create(logFile)
sink(logFile, append = T, split=T)

metaData <- read.table(metaDataFile, header=T, sep="\t", stringsAsFactors=F)
metaData <- metaData[,c("Sample.name","Population.code","Superpopulation.code")]

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

# Get the list of samples
sample.id <- getScanID(genoData)
thousandGenomID <- sample.id[!is.na(str_extract(sample.id, "^H|N"))]
Gen3GID <- sample.id[!(sample.id %in% thousandGenomID)]
cat(paste0("Number of thousand genomes samples : ", length(thousandGenomID), "\n"))
cat(paste0("Number of Gen3G samples : ", length(Gen3GID), "\n"))
#/---------------------------------------------------------------------------------------------------------------------#

#\------------------------------------------- KING-robust kinship estimation ------------------------------------------#
cat("########################################\nKING-robust kinship estimation\n########################################\n")
dir.create(file.path(savePath,"FirstIter"))

# Prune the SNPs. Using the suggested parameters.
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, ld.threshold=sqrt(0.1), maf=0.05, missing.rate=0.01, num.thread=10)
cat(paste0("Number of SNPs : ", length(unlist(snpset)), "\n"))

# Infer relatedness with the KING-robust method inplemented in SNPRelate
kingRobust <- snpgdsIBDKING(gds, snp.id=unlist(snpset, use.names=F), num.thread=10)

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

# Create annotation df
ancestry <- rbind(metaData, as.data.frame(list(Sample.name=Gen3GID, Population.code=rep("Gen3G",length(Gen3GID)), Superpopulation.code=rep("Gen3G",length(Gen3GID))),stringsAsFactors=F))
pcaAll <- merge(pcsFirstIter, ancestry, by.x=0, by.y=1)
pcaAll <- pcaAll[order(pcaAll$Superpopulation.code), ]

# Plot the PCs from the first iteration
ggplot(varProp, aes(x=Order,y=Var*100)) + geom_point() + labs(x="",y="Variance proportion (%)")
ggsave(file.path(savePath,"FirstIter","pcs_var_plot.jpg"))

ggplot(pcaAll,mapping=aes(V1,V2,color=Superpopulation.code)) + geom_point() + labs(x=paste0("PC1 (",round(varProp[1,1]*100, 3),"%)"),y=paste0("PC2 (",round(varProp[2,1]*100,3),"%)")) + scale_color_manual(values=rainbow(length(unique(pcaAll$Superpopulation.code))))
ggsave(file.path(savePath,"FirstIter","pcair_pc1_pc2.jpg"))

ggplot(pcaAll,mapping=aes(V1,V3,color=Superpopulation.code)) + geom_point() + labs(x=paste0("PC1 (",round(varProp[1,1]*100, 3),"%)"),y=paste0("PC3 (",round(varProp[3,1]*100,3),"%)")) + scale_color_manual(values=rainbow(length(unique(pcaAll$Superpopulation.code))))
ggsave(file.path(savePath,"FirstIter","pcair_pc1_pc3.jpg"))

ggplot(pcaAll,mapping=aes(V2,V3,color=Superpopulation.code)) + geom_point() + labs(x=paste0("PC2 (",round(varProp[2,1]*100, 3),"%)"),y=paste0("PC3 (",round(varProp[3,1]*100,3),"%)")) + scale_color_manual(values=rainbow(length(unique(pcaAll$Superpopulation.code))))
ggsave(file.path(savePath,"FirstIter","pcair_pc2_pc3.jpg"))

save(pcairFirstIter,pcsFirstIter,varProp,ancestry,file=PcairObjSavePath)
#/---------------------------------------------------------------------------------------------------------------------#

#\----------------------------------------------------- PC-Relate -----------------------------------------------------#
cat("########################################\nPC-Relate\n########################################\n")
PcrelateObjSavePath <- file.path(savePath,"SecondIter","pcrelate.RData")
dir.create(file.path(savePath,"SecondIter"))

# Create a GenomicDataIterator reading only the pruned SNPs
genoDataIter <- GenotypeBlockIterator(genoData, snpInclude=unlist(snpset, use.names=F))
# PC-Realte, Compute relatedness accounting for ancestry (first 4 PC of PC-AiR)
pcrel <- pcrelate(genoDataIter, pcs=pcairFirstIter$vectors[,1:4], training.set=pcairFirstIter$unrels)
# Extract a GRM, N x N kinship matrix to use in the second iteration of PC-AiR
pcrelGRM <- pcrelateToMatrix(pcrel)
# Plot the relatedness
ggplot(pcrel$kinBtwn,aes(x=1:length(ID1),y=kin,label=ifelse(kin>0.1,paste(ID1,ID2,sep=" - "),"")))+geom_point()+geom_text(hjust=0, vjust=0)+labs(x="",y="kinship")
ggsave(file.path(savePath,"SecondIter","kinship.jpg"))
ggplot(pcrel$kinBtwn,aes(x=k0,y=kin,label=ifelse(kin>0.1,paste(ID1,ID2,sep=" - "),"")))+geom_point()+geom_text(hjust=0, vjust=0)
ggsave(file.path(savePath,"SecondIter","kinship_IBD0.jpg"))
ggplot(pcrel$kinBtwn,aes(x=k0,y=k2,label=ifelse(k0<0.6,paste(ID1,ID2,sep=" - "),"")))+geom_point()+geom_text(hjust=0, vjust=0)
ggsave(file.path(savePath,"SecondIter","IBD0_IBD2.jpg"))
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
ggplot(varPropSecondIter, aes(x=Order,y=Var*100)) + geom_point() + labs(x="",y="Variance proportion (%)")
ggsave(file.path(savePath,"SecondIter","pcs_var_plot.jpg"))

# Create annotation df
pcaAllSecondIter <- merge(pcsSecondIter, ancestry, by.x=0, by.y=1)
pcaAllSecondIter <- pcaAllSecondIter[order(pcaAllSecondIter$Superpopulation.code),]

ggplot(pcaAllSecondIter,mapping=aes(V1,V2,color=Superpopulation.code)) + geom_point() + labs(x=paste0("PC1 (",round(varPropSecondIter[1,1]*100, 3),"%)"),y=paste0("PC2 (",round(varPropSecondIter[2,1]*100,3),"%)")) + scale_color_manual(values=rainbow(length(unique(pcaAllSecondIter$Superpopulation.code))))
ggsave(file.path(savePath,"SecondIter","pcair_pc1_pc2.jpg"))

ggplot(pcaAllSecondIter,mapping=aes(V1,V3,color=Superpopulation.code)) + geom_point() + labs(x=paste0("PC1 (",round(varPropSecondIter[1,1]*100, 3),"%)"),y=paste0("PC3 (",round(varPropSecondIter[3,1]*100,3),"%)")) + scale_color_manual(values=rainbow(length(unique(pcaAllSecondIter$Superpopulation.code))))
ggsave(file.path(savePath,"SecondIter","pcair_pc1_pc3.jpg"))

ggplot(pcaAllSecondIter,mapping=aes(V2,V3,color=Superpopulation.code)) + geom_point() + labs(x=paste0("PC2 (",round(varPropSecondIter[2,1]*100, 3),"%)"),y=paste0("PC3 (",round(varPropSecondIter[3,1]*100,3),"%)")) + scale_color_manual(values=rainbow(length(unique(pcaAllSecondIter$Superpopulation.code))))
ggsave(file.path(savePath,"SecondIter","pcair_pc2_pc3.jpg"))

save(pcairSecondIter,pcsSecondIter,varPropSecondIter,pcaAllSecondIter,file=PcairSecItObjSavePath)
#/---------------------------------------------------------------------------------------------------------------------#

#\--------------------------------------------------- Find outliers ---------------------------------------------------#
cat("########################################\nFinding outliers \n########################################\n")
cat("Using CEU and TSI populations as reference, and a threshold of 5 SD from the mean\n")
dir.create(file.path(savePath,"Outliers"))

# Build reference
pcaEuro <- subset(pcaAllSecondIter, Population.code %in% c("CEU","TSI"))
pcaEuroSdMean <- reduce(apply(pcaEuro[,2:6], 2, function(pc) as.data.frame(list(Mean=mean(pc),SD=sd(pc)))),rbind)

# Test Gen3G samples
gen3G_samples <- subset(pcaAllSecondIter, Population.code == "Gen3G", select=1:6)
gen3G_samples_dist <- sapply(1:5, function(pc) (gen3G_samples[, pc+1]-pcaEuroSdMean$Mean[pc])/pcaEuroSdMean$SD[pc])
rownames(gen3G_samples_dist) <- gen3G_samples$Row.names
colnames(gen3G_samples_dist) <- c("SD_PC1", "SD_PC2", "SD_PC3", "SD_PC4", "SD_PC5")

# Select outliers
gen3G_outliers <- gen3G_samples_dist[apply(gen3G_samples_dist, 1, function(x) any(abs(x) >= 5)), ]
write.table(gen3G_outliers, file=outliers_file)
write.table(gen3G_samples_dist, file=file.path(savePath, "Outliers", "SD_from_mean_all.tsv"))

# Plot PCA w/ outliers
pcaAllSecondIter[pcaAllSecondIter$Row.names %in% rownames(gen3G_outliers), ]$Superpopulation.code <- "Gen3G_outliers"
pcaAllSecondIter <- pcaAllSecondIter[order(pcaAllSecondIter$Superpopulation.code),]

ggplot(pcaAllSecondIter,mapping=aes(V1,V2,color=Superpopulation.code)) + geom_point() + labs(x=paste0("PC1 (",round(varPropSecondIter[1,1]*100, 3),"%)"),y=paste0("PC2 (",round(varPropSecondIter[2,1]*100,3),"%)")) + scale_color_manual(values=rainbow(length(unique(pcaAllSecondIter$Superpopulation.code))))
ggsave(file.path(savePath,"Outliers","pcair_pc1_pc2.jpg"))

ggplot(pcaAllSecondIter,mapping=aes(V1,V3,color=Superpopulation.code)) + geom_point() + labs(x=paste0("PC1 (",round(varPropSecondIter[1,1]*100, 3),"%)"),y=paste0("PC3 (",round(varPropSecondIter[3,1]*100,3),"%)")) + scale_color_manual(values=rainbow(length(unique(pcaAllSecondIter$Superpopulation.code))))
ggsave(file.path(savePath,"Outliers","pcair_pc1_pc3.jpg"))

ggplot(pcaAllSecondIter,mapping=aes(V2,V3,color=Superpopulation.code)) + geom_point() + labs(x=paste0("PC2 (",round(varPropSecondIter[2,1]*100, 3),"%)"),y=paste0("PC3 (",round(varPropSecondIter[3,1]*100,3),"%)")) + scale_color_manual(values=rainbow(length(unique(pcaAllSecondIter$Superpopulation.code))))
ggsave(file.path(savePath,"Outliers","pcair_pc2_pc3.jpg"))
#/---------------------------------------------------------------------------------------------------------------------#

#\---------------------------------------------------- Close files ----------------------------------------------------#
snpgdsClose(gds)
close(genoData)
sink()
#/---------------------------------------------------------------------------------------------------------------------#