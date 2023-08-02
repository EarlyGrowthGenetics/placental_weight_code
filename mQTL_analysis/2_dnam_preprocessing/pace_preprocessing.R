# Attach package
library(PACEanalysis)
library(reshape)

# Load data
args = commandArgs(trailingOnly=TRUE)
load(args[1]) # RGset_ext_combined
savePath <- args[2]
allphenodata<-read.table(args[3],header=TRUE,sep = ";",dec=",",as.is = T,na.strings = "")

# Exploratory data analysis
EDAresults<-ExploratoryDataAnalysis(RGset=RGset_ext_combined,
                                    globalvarexplore=c("BWT","Sex"),
                                    DetectionPvalMethod="SeSAMe",
                                    DetectionPvalCutoff=0.05,
                                    minNbeads=3,
                                    FilterZeroIntensities=TRUE,
                                    destinationfolder=savePath,
                                    savelog=TRUE,
                                    cohort="Gen3G",analysisdate="20210709")

save(EDAresults, file=file.path(savePath, "EDAresults2.RData"))

# Check whether reported sex matches inferred sex ()
rawbetas<-preprocessRaw(RGset_ext_combined)
rawbetas<-mapToGenome(rawbetas)
estSex <- getSex(rawbetas)
rawbetas <- addSex(rawbetas, sex = estSex)
png(file.path(savePath, "Gen3G_20210713_Signal_Intensities_on_X_and_Y_Chromosomes_by_Sex_CA.png"))
plotSex(rawbetas,id=pData(rawbetas)$ID)
dev.off()

# Identify bad samples
SexWrongRemove<-c()
TooManyFailedProbes<-c("203519500044_R04C01")
PossibleContamination<-c("200598360050_R07C01","203519500053_R05C01")
Unintentionalreps<-c()
allsamplestoexclude<-c(SexWrongRemove,TooManyFailedProbes,PossibleContamination,Unintentionalreps)


# Pre-processing
processedOut<-preprocessingofData(RGset=RGset_ext_combined,
                                  SamplestoRemove=allsamplestoexclude,
                                  ProbestoRemove=EDAresults$ProbestoRemove,
                                  destinationfolder=savePath,
                                  compositeCellType="Placenta",
                                  KchooseManual=NULL,
                                  savelog=TRUE,
                                  cohort="Gen3G",analysisdate="20210709")

write.table(processedOut$Omega, file = file.path(savePath, "estimated_cell_types.tsv"), sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

betasabovedetection<-detectionMask(processedBetas=processedOut$processedBetas,
                                   DetectionPvals=EDAresults$DetectionPval,
                                   DetectionPvalCutoff=0.05,
                                   IndicatorGoodIntensity=EDAresults$IndicatorGoodIntensity,
                                   destinationfolder=savePath,
                                   cohort="Gen3G",analysisdate="20210709")

Betasnooutliers<-outlierprocess(processedBetas=betasabovedetection,
                                quantilemethod="EmpiricalBeta",
                                trimming=FALSE,
                                pct=0.005,
                                destinationfolder=savePath,
                                cohort="Gen3G",analysisdate="20210709")
