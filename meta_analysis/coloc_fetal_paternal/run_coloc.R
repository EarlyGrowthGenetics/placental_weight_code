library(coloc)
library(tidyverse)
args<-commandArgs(TRUE)
i<-args[1]
analysis<-args[2]
if(analysis %in% c("sex")){
  data<-read_tsv(paste("locus_",i,"_joined",sep=""))
}else if(analysis %in% c("gest")){
  data<-read_tsv(paste("gest_locus_",i,"_joined",sep=""))
}
my.res<-coloc.abf(dataset1=list(beta=data$p1_beta,varbeta=data$p1_se^2,N=max(data$p1_n),sdY=1,MAF=data$eaf,type="quant"),dataset2=list(beta=data$p2_beta,varbeta=data$p2_se^2,N=max(data$p2_n),sdY=1,MAF=data$p2_eaf,type="quant"))
sink(paste("locus_",analysis,"_",i,"_out",sep=""))
print(my.res)
sink()
