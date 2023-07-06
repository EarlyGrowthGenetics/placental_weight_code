
library(data.table)
library(dplyr)
library(ggplot2)
library(ggtext)
library(tidyverse)
library(patchwork)

# Set the path to the folder containing the .txt files
path <- '/home/christopher/Desktop/child_gest/LDSC_Files/Fetal/cleaned_data'

# Get a list of all .txt files in the folder
txt_files <- list.files(path, pattern = "\\.txt$", full.names = TRUE)

# Loop through each .txt file, read it into a dataframe, and assign the file name as the dataframe name
for (file in txt_files) {
  df <- read.table(file, header = TRUE, sep = "\t")
  file_name <- gsub("\\.txt$", "", basename(file))
  assign(file_name, df)
}

df_list <- list(ALSPAC, DNBC, EFSOCH, FS, GENR, GOYA_control_mothers, GOYA_obese_children, GOYA_obese_mothers,
                GSAV2, IHPS_casesCtrls_MEGAex, IHPS_CIDR, INMAGSA, INMAOmni, iPSYCH, Meta, MoBa,
                NFBC1966, NFBC1986, OPI, PANIC, RaineStudy, Roskilde)

names(df_list) <- c("ALSPAC", "DNBC", "EFSOCH", "FS", "GENR", "GOYA_control_mothers", "GOYA_obese_children", "GOYA_obese_mothers", 
                    "GSAV2", "IHPS_casesCtrls_MEGAex", "IHPS_CIDR", "INMAGSA", "INMAOmni", "iPSYCH", "Meta", "MoBa", "NFBC1966", "NFBC1986", 
                    "OPI", "PANIC", "RaineStudy", "Roskilde")


# Load heterogeneity results from the meta

ht <- fread('/home/christopher/Desktop/child_gest/LDSC_Files/Fetal/cleaned_data/het.csv')

ht <- ht %>% rename("I^2" = Isq,
                    'Het P' = Het_P)


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

df$study <- recode(df$study, 'GOYA_control_mothers' = "GOYA_ctr_mum", 
                                 'GOYA_obese_children' = "GOYA_ob_child",
                                 'GOYA_obese_mothers' = "GOYA_ob_mum",
                                 'IHPS_casesCtrls_MEGAex' = "IHPS_MEGAex")

# create a list of dataframes for each group
groups <- split(df, df$SNP)


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
plot1 <- rs1021508 + rs10486660 + rs10925945 + rs112635299 + rs11708067 + rs11756568 
plot2 <- rs11866404 + rs12529634 + rs12543725 + rs138715366 + rs140691414 + rs1434836
plot3 <- rs150138294 + rs1655296 + rs1801253 + rs2237892 + rs3822394 + rs4953353
plot4 <- rs541641049 + rs55958435 + rs57790054 + rs6040436 + rs6078190 + rs6456014
plot5 <- rs6557677 + rs67265526 + rs7177338 + rs723177 + rs72801474 + rs74457440
plot6 <- rs7722058 + rs7783810 + rs876987 +rs9800506 + rs9817452 + plot_layout(ncol = 3)

plots <- list(plot1, plot2, plot3, plot4, plot5, plot6)

# save plots
for (i in 1:length(plots)) {
  ggsave(
    filename = paste0("/home/christopher/placental_weight_code/Heterogeneity/Fetal/Results/Meta_Analysis_Forest_Plots/top_hits_forest", i, ".pdf"),
    plot = plots[[i]],
    width = 14.5,
    height = 12,
    device = "pdf",
    dpi = 600
  )
}


for (i in 1:length(plots)) { 
ggsave(path = '/home/christopher/placental_weight_code/Heterogeneity/Fetal/Results/Meta_Analysis_Forest_Plots',
       filename = paste0("top_hits_forest", i, ".png"),
       plot = plots[[i]], 
       width = 14.5, 
       height = 12, 
       device = "png", 
       dpi = 900)
}

for (i in 1:length(plots)) {
  ggsave(
    filename = paste0("top_hits_forest", i, ".eps"),
    plot = plots[[i]],
    width = 14.5,
    height = 12,
    device = "eps",
    dpi = 600,
    path = '/home/christopher/placental_weight_code/Heterogeneity/Fetal/Results/Meta_Analysis_Forest_Plots'
  )
}


# df <- arrange(df, desc(SNP))
# 
# write.table(df, '/home/christopher/Desktop/child_gest/LDSC_Files/Fetal/cleaned_data/combined_data.txt', col.names=T, row.names=F, quote=F,sep='\t')
# 

