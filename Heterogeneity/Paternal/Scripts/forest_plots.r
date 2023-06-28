
library(data.table)
library(dplyr)
library(ggplot2)
library(ggtext)
library(tidyverse)
library(patchwork)

# Set the path to the folder containing the .txt files
path <- '/home/christopher/Desktop/child_gest/LDSC_Files/Paternal/Cleaned_data'

# Get a list of all .txt files in the folder
txt_files <- list.files(path, pattern = "\\.txt$", full.names = TRUE)

# Loop through each .txt file, read it into a dataframe, and assign the file name as the dataframe name
for (file in txt_files) {
  df <- read.table(file, header = TRUE, sep = "\t")
  file_name <- gsub("\\.txt$", "", basename(file))
  assign(file_name, df)
}

df_list <- list(DBDS, EFSOCH, GDAffy, GDIllum, GSAV2, HUNT, Meta, MoBa, RaineStudy)

names(df_list) <- c("DBDS", "EFSOCH", "GDAffy", "GDIllum", "GSAV2", "HUNT", "Meta", "MoBa", "RaineStudy")

# extract a list of SNPS to use to get the heterogeneity stats
SNP <- Meta %>% select(SNP) 

# Load meta results and extract heterogeneity results from the meta

ht <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/Paternal/Clean/pw_paternal_sex_gest.gz')

ht <- ht %>% rename(SNP = MarkerName)

ht <- inner_join(SNP, ht) %>% select(SNP, HetISq, HetPVal)

ht <- ht %>% rename("I^2" = HetISq,
                    'Het P' = HetPVal)


# Iterate over the list of dataframes and create column called study to identify each study
for (i in 1:length(df_list)) {
  df_list[[i]]$study <- names(df_list)[i]
}


# combine all dataframes into one
df <- bind_rows(df_list)

df$ci.lb <- df$BETA - 1.96*df$SE    # calculate lower bounds of 95% CI
df$ci.ub <- df$BETA + 1.96*df$SE    # calculate upper bounds of 95% CI

# merge the heterogeneity results into dataframe

df <- inner_join(df, ht, by='SNP')

# Make p values 3 decimal places

df$'Het P' <- format(round(df$'Het P', 3), nsmall = 3)


# create a list of dataframes for each group
groups <- split(df, df$SNP)


# make the labels for the meta red
a <- ifelse(df$study == "Meta", "red", "black")



# loop through each group and create a forest plot after ordering the studies on the y axis by sample size

for (g in names(groups)) {
  plot <- ggplot(data = groups[[g]], aes(x = BETA, y = fct_reorder(study, N, .desc = TRUE), xmin = ci.lb, xmax = ci.ub, color = study)) +
    geom_pointrange() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0) +
    ggtitle(paste0("", g, 
                   paste0('    I2 = ', {groups[[g]]$'I^2'}), 
                   paste0(",   HET P = ", {groups[[g]]$'Het P'}))) +
    xlab("Beta (95% CI)") +
    ylab("Study") +
    theme(axis.title.y = element_text(size = 12)) +
    scale_color_manual(values = c("Meta" = "red")) +
    scale_y_discrete(labels = function(x) ifelse(x == "Meta", paste0("<span style='color:red'>", x, "</span>"), x)) +
    theme(axis.text.y = element_markdown(size = 10, face = "bold")) +
    theme(plot.title = element_text(size = 10, face = "bold")) +
    guides(color = "none") +
    theme(panel.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  assign(g, plot, envir = .GlobalEnv)
}

# combine plots
plot1 <- rs2207099 + rs75512885 

# save plots

ggsave(
  filename = "paternal_top_hits_forest.png",
  plot = plot1, 
  width = 14.5, 
  height = 12, 
  dpi = 600,
  path = "/home/christopher/placental_weight_code/Heterogeneity/Paternal/Results/Meta_Analysis_Forest_Plots",
  device = "png"
)



# df <- arrange(df, desc(SNP))
# 
# write.table(df, '/home/christopher/Desktop/child_gest/LDSC_Files/Paternal/Cleaned_data/combined_data.txt', col.names=T, row.names=F, quote=F,sep='\t')


