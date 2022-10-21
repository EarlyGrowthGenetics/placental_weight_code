#! /usr/bin/env Rscript
#
# Jess Tyrrell
# Genetics of Complex Traits
# University of Exeter Medical School
# Email: J.Tyrrell@exeter.ac.uk
# Date:  16/03/2017
# R Script to run various MR analyses using individual SNPs for BOLT and BGENIE


# Define input arguments
args <- commandArgs(TRUE)

test_trait <- args[1] #name of trait
snp_info <- args[2] #location of snp_info_file
summ_stats_file <- args[3] #location of file containing individual BetaYG and seBetaYG for a range of traits
output_dir <- args[4] #directory in which to save results and plots

if(identical(substr(output_dir,nchar(output_dir),nchar(output_dir)),"/")) {
   output_dir <- substr(output_dir,1,nchar(output_dir)-1)
}


# Check all required arguments are present - if not display usage information
if(length(args)<4) {
   system(paste("echo -e \"\nIncorrect number of arguments given (",length(args)," given, 4 expected). Correct usage:\"",sep=""))
   system("echo -e \"Rscript general_iv_script.R [TRAIT NAME] [SNP INFO FILE] [SUMMARY STATS FILE] [OUTPUT DIRECTORY]\n\"")
   q("no", status = 1, runLast = FALSE)
}

if(!file.exists(snp_info)) {
   system(paste("echo -e \"\nSNP info file \\\"",snp_info,"\\\" does not exist... exiting!\n\"",sep=""))
   q("no", status = 1, runLast = FALSE)
}
if(!file.exists(summ_stats_file)) {
   system(paste("echo -e \"\nSummary stats file \\\"",summ_stats_file,"\\\" does not exist... exiting!\n\"",sep=""))
   q("no", status = 1, runLast = FALSE)
}

if(!file.exists(output_dir)) {
   system(paste("echo -e \"\nOutput directory \\\"",output_dir,"\\\" does not exist... exiting!\n\"",sep=""))
   q("no", status = 1, runLast = FALSE)
}





#Load relevant packages
library(Cairo)
library(plyr)


#Load data into R
snp_info<-read.table(snp_info, h=TRUE, sep="\t", strip.white=TRUE)
snp_info_specific<-snp_info[which(snp_info$Test_trait==test_trait),]
results<-read.table(summ_stats_file, h=TRUE, sep="\t", strip.white=TRUE)

#Merge the snp_info and results file and flip to trait raising allele
stats <-merge(results, snp_info_specific, by="SNP")
stats$BetaYG<-ifelse(as.character(stats$Trait_raising)==as.character(stats$A1), stats$BetaYG, -stats$BetaYG)

# Run MR 
# Multiple outcomes of interest - use code below (will create table of results that can be used for plots)
# Generate new variables
stats$betaIV       = stats$BetaYG/stats$BetaXG
stats$weights      = (stats$seBetaYG/stats$BetaXG)^-2
# Run IV analysis across the traits
betaIVW<-ddply(stats, "Trait", transform, betaIVW=(sum(BetaYG*BetaXG*seBetaYG^-2)/sum(BetaXG^2*seBetaYG^-2)))

betaIVW2<-ddply(stats, "Trait", summarize, betaIVW2=summary(lm(BetaYG~BetaXG-1, weights=seBetaYG^-2))$coeff[1,1])
sebetaIVW2 <-ddply(stats, "Trait", summarize, sebetaIVW2=summary(lm(BetaYG~BetaXG-1, weights=seBetaYG^-2))$coeff[1,2]/
					min(summary(lm(BetaYG~BetaXG-1, weights=seBetaYG^-2))$sigma, 1))
# Run Egger analysis across the traits
betaEgger<-ddply(stats, "Trait", summarize, betaEgger=summary(lm(BetaYG~BetaXG, weights=seBetaYG^-2))$coeff[2,1])
sebetaEgger <-ddply(stats, "Trait", summarize, sebetaEgger=summary(lm(BetaYG~BetaXG, weights=seBetaYG^-2))$coeff[2,2]/
					min(summary(lm(BetaYG~BetaXG, weights=seBetaYG^-2))$sigma, 1))
egger_int<-	ddply(stats, "Trait", summarize, egger_int=summary(lm(BetaYG~BetaXG, weights=seBetaYG^-2))$coeff[1,1])
int_p <-	ddply(stats, "Trait", summarize, int_p=summary(lm(BetaYG~BetaXG, weights=seBetaYG^-2))$coeff[1,4])										
# Combine IV stat with main data for the median analyses
stats<-cbind(stats, betaIVW$betaIVW) 					
stats$penalty      = pchisq(stats$weights*(stats$betaIV-stats$betaIVW)^2, df=1, lower.tail=FALSE)
stats$pen.weights  = stats$weights*pmin(1, stats$penalty*20)               # penalized weights
# Functions for the median analyses
weighted.median <- function(betaIV.in, weights.in) {
  betaIV.order  = betaIV.in[order(betaIV.in)]
  weights.order = weights.in[order(betaIV.in)]
weights.sum = cumsum(weights.order)-0.5*weights.order
weights.sum = weights.sum/sum(weights.order)
      below = max(which(weights.sum<0.5))
  weighted.est  = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
                  (0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
  return(weighted.est) }
  
weighted.median.boot = function(BetaXG.in, betaYG.in, seBetaXG.in, sebetaYG.in, weights.in){
med = NULL
for(i in 1:1000){
 BetaXG.boot = rnorm(length(BetaXG.in), mean=BetaXG.in, sd=seBetaXG.in)
 betaYG.boot = rnorm(length(betaYG.in), mean=betaYG.in, sd=sebetaYG.in)
 betaIV.boot  = betaYG.boot/BetaXG.boot
 med[i] = weighted.median(betaIV.boot, weights.in)
 }
return(sd(med)) }
# Run Median IV analyses
betaWM<-ddply(stats, "Trait", summarize, betaWM=(weighted.median(betaIV, weights)))
sebetaWM<-ddply(stats, "Trait", summarize, sebetaWM=(weighted.median.boot(BetaXG, BetaYG, seBetaXG, seBetaYG, weights)))
betaPWM<-ddply(stats, "Trait", summarize, betaPWM=(weighted.median(betaIV, pen.weights)))
sebetaPWM<-ddply(stats, "Trait", summarize, sebetaPWM=(weighted.median.boot(BetaXG, BetaYG, seBetaXG, seBetaYG, pen.weights)))

# Collate data into one table
MR_results<-cbind(betaWM, sebetaWM, betaPWM, sebetaPWM, betaIVW2, sebetaIVW2, betaEgger, sebetaEgger, egger_int, int_p)
MR_results_final<-subset(MR_results, select = c("Trait", "betaWM", "sebetaWM", "betaPWM", "sebetaPWM", "betaIVW2", "sebetaIVW2", "betaEgger", "sebetaEgger", "egger_int", "int_p"))
MR_results_final$tWM= abs(MR_results_final$betaWM/MR_results_final$sebetaWM)                                                       
MR_results_final$tPWM= abs(MR_results_final$betaPWM/MR_results_final$sebetaPWM)
MR_results_final$tIVW= abs(MR_results_final$betaIVW2/MR_results_final$sebetaIVW2)
MR_results_final$tegger=abs(MR_results_final$betaEgger/MR_results_final$sebetaEgger)			
# Calculate p values - t-distribution for IVW and Egger, normal distribution for median analyses
n_snp = length(unique(stats$SNP))
MR_results_final$pIVW=2*pt((-MR_results_final$tIVW), df=n_snp-1)
MR_results_final$pEgger=2*pt((-MR_results_final$tegger), df=n_snp-2)
MR_results_final$pWM=2*pnorm((-MR_results_final$tWM))
MR_results_final$pPWM=2*pnorm((-MR_results_final$tPWM))
MR_results_final$IVR_int = 0
MR_results_final$n_snp = n_snp
MR_results_final<-MR_results_final[c("Trait", "betaIVW2", "sebetaIVW2", "tIVW", "pIVW", "IVR_int", "betaEgger", "sebetaEgger", "tegger", "pEgger", "egger_int", "int_p", "betaWM", "sebetaWM", "tWM", "pWM", "betaPWM", "sebetaPWM", "tPWM", "pPWM", "n_snp" )]
#Output results
output_tablename<-paste(output_dir,"/IV_results_",test_trait,".txt",sep="")
write.table(MR_results_final,output_tablename,sep="\t",row.names=FALSE)

#Plot results
# Load in SNP association data and summary stats

data <- stats
rm(stats)
stats <- MR_results_final

# Loop through traits
for (traitname in unique(unlist(data$Trait))) {
cur_data  <- subset(data,Trait==traitname)
cur_stats <- subset(stats,Trait==traitname)


# set strings for filename and titles
output_filename <- paste(output_dir,"/",test_trait,"_SNP_",traitname,".png",sep="")
ylabel <- paste(traitname,test_trait," SNP Beta")
xlabel <- paste(test_trait," Beta")
title  <- paste(test_trait," SNPs vs. ",traitname," Betas")


# open plotting file
CairoPNG(output_filename,800,800)


plot(cur_data$BetaXG,cur_data$BetaYG,xlim=range(c(min(0,cur_data$BetaXG-1.96*cur_data$seBetaXG),max(0,cur_data$BetaXG+1.96*cur_data$seBetaXG))),ylim=range(c(min(0,cur_data$BetaYG-1.96*cur_data$seBetaYG),max(0,cur_data$BetaYG+1.96*cur_data$seBetaYG))),pch=19,xlab=xlabel,ylab=ylabel,main=title)

arrows(cur_data$BetaXG,cur_data$BetaYG-1.96*cur_data$seBetaYG, cur_data$BetaXG, cur_data$BetaYG+1.96*cur_data$seBetaYG, length=0.05, angle=90, code=3)

arrows(cur_data$BetaXG-1.96*cur_data$seBetaXG,cur_data$BetaYG, cur_data$BetaXG+1.96*cur_data$seBetaXG, cur_data$BetaYG, length=0.05, angle=90, code=3)

abline(h=0)
abline(v=0)

abline(cur_stats$egger_int,cur_stats$betaEgger,col="red")
abline(cur_stats$IVR_int,cur_stats$betaIVW2,col="blue")
abline(cur_stats$IVR_int, cur_stats$betaWM, col="deeppink")
abline(cur_stats$IVR_int, cur_stats$betaPWM, col="green4")

legend("topright",c(paste("MR-Egger: P = ",format(cur_stats$pEgger, scientific = TRUE, digits = 2)),paste("IVW: P = ",format(cur_stats$pIVW, scientific = TRUE, digits = 2)),paste("Median IV: P = ",format(cur_stats$pWM, scientific = TRUE, digits = 2)),paste("Penalised Median IV: P = ",format(cur_stats$pPWM, scientific = TRUE, digits = 2))),col=c("red","blue","deeppink","green4"),lty=c(1,1,1,1))


dev.off()

}


