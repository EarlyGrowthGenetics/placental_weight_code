#\---------------------------------------------------- Description ----------------------------------------------------#
# Filters out CpG with low variance
#/---------------------------------------------------------------------------------------------------------------------#

#\----------------------------------------------------- Libraries -----------------------------------------------------#
library(ggplot2)
library(gridExtra)
#/---------------------------------------------------------------------------------------------------------------------#

#\----------------------------------------------------- Load data -----------------------------------------------------#
args = commandArgs(trailingOnly=TRUE)

load(args[1]) # betas_filtered
varThld <- args[2]
logFile <- args[3]
savePath <- args[4]
savePath_samples <- args[5]
plots_savePath <- args[6]

file.create(logFile)
sink(logFile, append = T, split=T)
#/---------------------------------------------------------------------------------------------------------------------#

#\----------------------------------------------- Plot var distribution -----------------------------------------------#
cat("########################################\nCalculate variance\n########################################\n")
mean_var <- as.data.frame(list(mean=apply(betas_filtered,1,function(x) mean(x, na.rm=T)),var=apply(betas_filtered,1,function(x) var(x, na.rm=T))))
summary(mean_var$var)
quantile(mean_var$var, c(0.01,0.05,0.1))

var_hist <- ggplot(mean_var, aes(x=var)) + geom_histogram()
var_hist_cut <- ggplot(mean_var, aes(x=var)) + geom_histogram() + xlim(c(0,0.0005))
var_density <- ggplot(mean_var, aes(x=mean, y=var)) + geom_bin2d(bins=100) + scale_fill_continuous(type = "viridis")
plots <- grid.arrange(grobs=list(var_hist,var_hist_cut,var_density), width=c(1,1), layout_matrix=rbind(c(1,2),c(3,3),c(3,3)))
ggsave(plots, file=plots_savePath)
#/---------------------------------------------------------------------------------------------------------------------#

#\------------------------------------------------------- Filter ------------------------------------------------------#
cat("########################################\nFiltering data\n########################################\n")
print(paste("Using variance threshold :", varThld))
betas_filtered_var <- subset(betas_filtered, mean_var$var >= varThld)
dim(betas_filtered_var)
print(paste0("The filter removed ", round((1-nrow(betas_filtered_var)/nrow(betas_filtered))*100, 2), "% of the probes"))
#/---------------------------------------------------------------------------------------------------------------------#

#\-------------------------------------------------------- Save -------------------------------------------------------#
save(betas_filtered_var, file=savePath)
write.table(colnames(betas_filtered_var), file = savePath_samples, sep="\t", row.names=F, col.names=F, quote = F)
sink()
#/---------------------------------------------------------------------------------------------------------------------#