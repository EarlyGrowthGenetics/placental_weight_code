
library(data.table)
library(dplyr)
library(metafor)
library(ggplot2)

# Set the path to the folder containing the .txt files
path <- '/home/christopher/Desktop/child_gest/LDSC_Files/cleaned_data'

# Get a list of all .txt files in the folder
txt_files <- list.files(path, pattern = "\\.txt$", full.names = TRUE)

# Loop through each .txt file, read it into a dataframe, and assign the file name as the dataframe name
for (file in txt_files) {
  df <- read.table(file, header = TRUE, sep = "\t")
  file_name <- gsub("\\.txt$", "", basename(file))
  assign(file_name, df)
}

# remove Meta analysis from the combined data as not wanted as part of the random effects meta

combined_data <- combined_data %>% filter(study!="Meta")



# Initialize empty dataframe to store results
results_df <- data.frame(SNP = character(), Beta = numeric(), SE = numeric(), P =numeric())


# Loop over variable x
for (i in combined_data$SNP) {
  # Subset data for current value of SNP
  dat_subset <- subset(combined_data, SNP == i)
  
  # Sort data frame based on sample size
  dat_subset <- dat_subset[order(dat_subset$N, decreasing = FALSE),]
  
  # Fit random effects model
  res <- rma(yi=BETA, vi=SE,  data=dat_subset, method="REML")
  
  # Plot forest plot
  forest_plot <- forest(res, slab = paste(unique(dat_subset$study)))
  
  # Save forest plot as png file
  pdf(file=paste0("/home/christopher/Desktop/child_gest/LDSC_Files/cleaned_data/Random_effects/Forest_plots/forest_plot_", i, ".pdf"))
  forest(res, slab = paste(unique(dat_subset$study)))
  dev.off() 
  
  # Plot funnel plot
  funnel_plot <- funnel(res)
 
  # Save funnel plot as png file
  pdf(file=paste0("/home/christopher/Desktop/child_gest/LDSC_Files/cleaned_data/Random_effects/Funnel_plots/funnel_plot_", i, ".pdf"))
  funnel(res)
  dev.off()

  # define which SNP
  SNP <- i
  
  #Extract Estimate
  Beta <- summary(res)$beta
    
  # Extract SE
  SE <- summary(res)$se
    
  # Extract P Value
  P <- summary(res)$pval
    
  # # Extract tau2
  tau2 <- summary(res)$tau2

  # Extract I^2 (inconsistency index)
  I2 <- summary(res)$I2
  
  # Extract H^2 (heterogeneity variance)
  H2 <- summary(res)$H2
  
  # Extract Q (Cochran's Q statistic)
  Q <- summary(res)$QE
  
  #Extract QP
  QP <- summary(res)$QEp
  
  # Create a new dataframe with the extracted values and current SNP
  temp_df <- data.frame(SNP = SNP, Beta = Beta, SE = SE, P = P, tau2 = tau2, I2 = I2, H2 = H2, Q = Q, Q_Pval = QP)

  # Append the new dataframe to the results_df
  results_df <- unique(rbind(results_df, temp_df))

}

results_df <-  results_df %>% select(SNP, Beta, SE, P, Q, Q_Pval)

# Define the columns that you want to round
cols_to_round <- c("Beta", "SE", "P", "Q", "Q_Pval")

# Round specific columns in the data frame to 3 decimal places
results_df[,cols_to_round] <- lapply(results_df[,cols_to_round], round, 3)

# extract fixed effects results

fe <-  fread('/home/christopher/Desktop/child_gest/LDSC_Files/cleaned_data/Meta.txt')

fe <-  fe %>% select(SNP, BETA, SE, P, N) %>% rename(FE_BETA = BETA, FE_SE = SE, FE_P = P)

results_df <- results_df %>% select(SNP, Beta, SE, P, Q, Q_Pval) %>% rename(RE_BETA = Beta, RE_SE = SE, RE_P = P)

results_df <- inner_join(fe, results_df, by='SNP')

results_df <- results_df %>% select("SNP", "FE_BETA", "FE_SE", "FE_P", "RE_BETA", "RE_SE", "RE_P", "Q", "Q_Pval", "N")
  
write.table(results_df, '/home/christopher/Desktop/child_gest/LDSC_Files/cleaned_data/Random_effects/Random_Effects_Results.txt',
            col.names=T, row.names=F, quote=F, sep='\t')

