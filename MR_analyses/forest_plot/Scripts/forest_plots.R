library(here)
library(glue)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)

df <- read.table(
  file = file.path(here(), 'MR_analyses/forest_plot/Data/plots_data.csv.gz'),
  header = T,
  sep = ","
  )


############## Plots

df$person <- paste(df$Genome, df$Analysis)

df$exposure[df$exposure == "Insulin Secretion"] <- "Disposition Index (Beta-cell Function)"
df$exposure[df$exposure == "Insulin Resistance"] <- "Fasting Insulin (Insulin Resistance)"


df$exposure <- factor(df$exposure, levels=c("Height","Fasting Glucose","Disposition Index (Beta-cell Function)",
                                            "Fasting Insulin (Insulin Resistance)", "SBP","DBP"))

#### Order data for plots

df <- df %>% arrange(exposure, Genome, Analysis)

df$person <- factor(df$person, levels = rev(unique(df$person)))

#### Convert BP to 10mm Hg instead of 1 mm Hg

df$Beta <- ifelse((df$exposure=="SBP" | df$exposure== "DBP"), df$Beta*10, df$Beta)
df$Se <- ifelse((df$exposure=="SBP" | df$exposure== "DBP"), df$Se*10, df$Se)

#### Change column names for legend

df$Analysis[df$Analysis=="Meta"] <- "Meta-analysis"
df$Analysis[df$Analysis=="WLM"] <- "WLM-adjusted"

###### Plot data


p1 <- ggplot(df, aes(Beta, person)) 

p1 <- p1 + aes(x=Beta, xmin=Beta-1.96*Se, xmax=Beta+1.96*Se, y=person, color=Genome) + geom_pointrange(aes(shape=Analysis), size = 1.5)  +
                xlim(-0.45,0.55) 
            
p1 <- p1 + facet_grid(exposure~.,scales="free",space="free",labeller=label_wrap_gen(width=2,multi_line=+T)) 

p1 <- p1 + geom_vline(xintercept=0) # add 0 line

p1 <- p1 + xlab("") + ylab("") + theme_bw() +
        theme(text = element_text(face = "bold", size=23, color = "black")) +
        scale_color_manual(values=c("#009E73","#D55E00")) + 
        xlab("Effect Size") + 
        theme(strip.text.y.right = element_text(angle=270)) +
        theme(axis.text.y.left = element_blank()) + 
        theme(axis.ticks.y.left = element_blank()) +
        theme(axis.text.x=element_text(size=20, face = "bold")) + 
        theme(legend.key = element_rect(fill = "white", colour = "black")) +
        guides(color = guide_legend(override.aes=list(fill = NA, linetype = 0, size=1), order=1),
               shape = guide_legend(override.aes = list(linetype = 0), order=2)) + 
        theme(legend.key.size = unit(1.5, 'cm')) +
        theme(legend.text = element_text(size=20)) 

         
  

# save the plot

figure_path <- file.path(here(), 'MR_analyses/forest_plot/Fig6_Mendelian_Randomization.png')
png(figure_path,
    width=1300, height=700)
p1
device <- dev.off()

