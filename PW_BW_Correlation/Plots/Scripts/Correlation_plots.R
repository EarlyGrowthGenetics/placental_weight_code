
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

df <- fread('BWT_Correlations/results/Combined_Results/combined_results.csv.gz')

df <- mutate (df,
      low = rg - (1.96*se), 
      high= rg + (1.96*se),
      h_low = h2_obs - (1.96*h2_obs_se), 
      h_high= h2_obs + (1.96*h2_obs_se),
    )


#### prepare trait1 for splitting into genome and analysis - 
#### levels will be used in plot legends

df$trait1 <- recode(df$trait1,
                    "fetal_ga_pwt" = "Fetal Meta",
                    "mum_ga" = "Maternal Meta",
                    "paternal_ga" = "Paternal Meta",
                    "wlm_child" = "Fetal WLM",
                    "wlm_maternal" = "Maternal WLM",
                    "wlm_paternal" = "Paternal WLM")


#### split trait1 into genome and analysis to be used in plot legends

df <- df %>% separate(trait1, c('Genome', 'Analysis'))


df$Analysis[df$Analysis=="Meta"] <- "Meta-analysis"
df$Analysis[df$Analysis=="WLM"] <- "WLM-adjusted"


#### Filter out pwt analyses

df1 <- filter(df, trait2!="Placenta Weight")

# keep only meta pw results for meta results for bw and wlm results together

df1 <- df1 %>% filter(Analysis=="Meta-analysis" & (trait2=="Offspring birth weight"|trait2=="Own birth weight") |
                                             Analysis=="WLM-adjusted" & (trait2=="BWT fetal effect"|trait2=="BWT maternal effect"))

#### create numeric variable to help order dataframe for plots

df1$or <- recode(df1$trait2,
                 "Own birth weight" = 1,
                 "Offspring birth weight" = 2,
                 "BWT fetal effect" = 3,
                 "BWT maternal effect" = 4)


#### order dataframe by trait2, then genome then analysis

df1 <- df1[
  with(df1, order(or, Genome, Analysis)),
]


df1$or <- NULL

#### Create numeric vector for correct order in plots

df1$n <- seq(1:12)


#### create a vector for x axis 

names <- c("Own Birth Weight \n (Fetal Genome)", "Offspring Birth Weight \n (Maternal Genome)", "Birth Weight \n Fetal WLM", "Birth Weight \n Maternal WLM")
names1 <- c("a", "b", "c", "d")

vlines <- c(3.5,6.5,9.5)
 

#### Plot genetic correlations

fp <- ggplot(data=df1, 
             aes(x=n, y=rg)) +
  scale_y_continuous(breaks = seq(-0.25, 1, by = 0.25)) +
  geom_rect(aes(xmin = (0+0.25),
                xmax = (3+0.5),
                ymin = -Inf, ymax = Inf, fill = '2002'), alpha = .2) + 
  geom_rect(aes(xmin = (4-0.5),
                xmax = (6+0.5),
                ymin = -Inf, ymax = Inf, fill = '2001'), alpha = .2) +
  geom_rect(aes(xmin = (7-0.5),
                xmax = (9+0.5),
                ymin = -Inf, ymax = Inf, fill = '2002'), alpha = .2) +
  geom_rect(aes(xmin = (10-0.5),
                xmax = (12+0.75),
                ymin = -Inf, ymax = Inf, fill = '2001'), alpha = .2) +
  geom_text(aes(x=c(0.5,3.75,6.75,9.75), label=names1), y=0.9, size=10, data=data.frame(names1)) + 
  geom_pointrange(aes(col=Genome, shape=Analysis, ymin=high, ymax=low, group=Genome), size=1) + 
  scale_color_manual(values = c("Fetal" = "#009E73",  
                                "Maternal" = "#D55E00", 
                                "Paternal" = "#0072B2"))  +
  scale_fill_manual(values=c("grey85","cornsilk1")) +
  geom_hline(yintercept=0, lty=2) + 
  geom_vline(xintercept=vlines, lty=2) +
  xlab("") + 
  ylab("Rg with PW (95% C.I.)")  +
  theme_bw() 

df1$n <- as.factor(df1$n)

fp1 <- fp + scale_x_discrete(limits=df1$n, breaks=c(2,5,8,11), 
                             label=as.character(names)) 
                             


fp1 <- fp1 + theme(plot.title = 
                      element_text(family = "Helvetica", face = "bold", size = (28)), 
                      axis.text.x = element_text(face = "bold", color="black", size = 14), 
                      axis.ticks.x = element_blank(),
                      axis.text.y = element_text(face = "bold", color = "black", size = 14),
                      axis.title.y = element_text(face = "bold", size = 14))

fp2 <- fp1 + guides(fill="none") 

fp2 <- fp2 + theme(legend.text=element_text(size=16), 
                   legend.title=element_text(size=16)) +
             theme(axis.title.y = element_text(size = 16, face="bold")) 

png(file="PW_BW_Correlation/Plots/Fig_4_pwt_bwt_correlation.png",
    width=900, height=500)
fp2
dev.off()



