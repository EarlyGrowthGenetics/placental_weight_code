library(coloc)
library(tidyverse)
library(magrittr)

# Read in arguments
args<-commandArgs(TRUE)
i<-args[1]

# read in the data
data<-read_tsv(paste("locus_",i,"_joined",sep=""))
data %<>% mutate(bw_se=abs(bw_beta/qnorm(bw_p/2))) %>% filter(bw_beta!=0)

# run coloc using original routine, allowing only 1 signal for each trait
my.res<-coloc.abf(dataset1=list(beta=data$pw_beta,varbeta=data$pw_se^2,N=max(data$pw_n),sdY=1,MAF=data$eaf,type="quant"),dataset2=list(beta=data$bw_beta,varbeta=data$bw_se^2,N=max(data$bw_n),sdY=1,MAF=data$bw_eaf,type="quant"))
## run coloc using susie routine to allow multiple signals
## NOTE: this requires LD, which we don't know for deCODE data so probably can't do this!!
## run susie on each trait
#pw<-susie()
## run coloc.susie on this output

# write out the data
sink(paste("locus_",i,"_out",sep=""))
print(my.res)
sink()
