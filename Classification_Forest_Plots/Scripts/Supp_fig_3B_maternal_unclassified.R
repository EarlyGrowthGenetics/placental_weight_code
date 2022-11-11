

library(dplyr)
library(data.table)
library(tidyverse)
library(patchwork)
library(ggpubr)

meta <- fread('meta/Classification_Forest_plots/Data/meta_ga_results_pos_beta.csv')

metaf <- subset(meta, select=c("chr", "pos", "rsid", "beta_fetal", "se_fetal"))
metam <- subset(meta, select=c("chr", "pos", "rsid", "beta_maternal", "se_maternal"))
metap <- subset(meta, select=c("chr", "pos", "rsid", "beta_paternal", "se_paternal"))

metaf <- metaf %>% rename( 
  Beta = beta_fetal,
  se = se_fetal,
)

metam <- metam %>% rename(
  Beta = beta_maternal,
  se = se_maternal,
)


metap <- metap %>% rename(
  Beta = beta_paternal,
  se = se_paternal,
)


metaf$Person <- paste("Fetal")
metam$Person <- paste("Maternal")
metap$Person <- paste("Paternal")


###### Load WLM data

wlm <- fread('meta/Classification_Forest_plots/Data/wlm_ga_results_pos_beta.csv', sep=",", head=TRUE)


wlmf <- subset(wlm, select=c("chr", "pos", "rsid", "wlm_beta_fetal", "wlm_se_fetal"))
wlmm <- subset(wlm, select=c("chr", "pos", "rsid", "wlm_beta_maternal", "wlm_se_maternal"))
wlmp <- subset(wlm, select=c("chr", "pos", "rsid", "wlm_beta_paternal", "wlm_se_paternal"))


wlmf <- wlmf %>% rename(
  Beta = wlm_beta_fetal,
  se = wlm_se_fetal,
)

wlmm <- wlmm %>% rename(
  Beta = wlm_beta_maternal,
  se = wlm_se_maternal,
)

wlmp <- wlmp %>% rename(
  Beta = wlm_beta_paternal,
  se = wlm_se_paternal,
)

wlmf$Person <- paste("Fetal WLM")
wlmm$Person <- paste("Maternal WLM")
wlmp$Person <- paste("Paternal WLM")


### Combine all data frames into one

df <- rbind(metaf,metam,metap,wlmf,wlmm,wlmp)


df <- filter(df, rsid=="rs150138294"|rsid=="rs10925945"|rsid=="rs4953353"|rsid=="rs11708067"|rsid=="rs74457440"|rsid=="rs72804545"|rsid=="rs75512885"|
                 rsid=="rs9800506"|rsid=="rs138715366"|rsid=="rs12543725"|rsid=="rs2237892"|rsid=="rs2168101"|rsid=="rs180435"|rsid=="rs7177338"|
                 rsid=="rs303998")



#### Create column for analysis type in plots

df$Analysis <- ifelse(grepl('WLM', df$Person), "WLM", "Meta")

df$Analysis[df$Analysis=="Meta"] <- "Meta-analysis"
df$Analysis[df$Analysis=="WLM"] <- "WLM-adjusted"

#### Create column for facet in plots

df$Genome[(grepl('Fetal', df$Person))] <- "Fetal"
df$Genome[(grepl('Maternal', df$Person))] <- "Maternal"
df$Genome[(grepl('Paternal', df$Person))] <- "Paternal"

##### Create gene column for facet labels

df$gene[df$rsid=="rs150138294"] <-"DCST2"
df$gene[df$rsid=="rs10925945"] <-"CHRM3"
df$gene[df$rsid=="rs4953353"] <-"EPAS1"
df$gene[df$rsid=="rs11708067"] <-"ADCY5"
df$gene[df$rsid=="rs74457440"] <-"PDLIM5"
df$gene[df$rsid=="rs72804545"] <-"EBF1"
df$gene[df$rsid=="rs75512885"] <-"EBF1"
df$gene[df$rsid=="rs9800506"] <-"FKBP5"
df$gene[df$rsid=="rs138715366"] <-"YKT6"
df$gene[df$rsid=="rs12543725"] <-"SLC45A4"
df$gene[df$rsid=="rs2237892"] <-"KCNQ1"
df$gene[df$rsid=="rs2168101"] <-"LMO1"
df$gene[df$rsid=="rs180435"] <-"SLC38A4"
df$gene[df$rsid=="rs7177338"] <-"FES"
df$gene[df$rsid=="rs303998"] <-"NLRP13"

#### Generate classification for background 

df$Classification[df$rsid=="rs72804545"|df$rsid=="rs2168101"|df$rsid=="rs180435"|df$rsid=="rs303998"] <- "Maternal"
df$Classification[df$rsid=="rs11708067"|df$rsid=="rs138715366"|df$rsid=="rs2237892"] <- "Fetal & Maternal"
df$Classification[df$rsid=="rs150138294"|df$rsid=="rs10925945"|df$rsid=="rs4953353"|df$rsid=="rs74457440"|df$rsid=="rs75512885"|
                  df$rsid=="rs9800506"|df$rsid=="rs12543725"|df$rsid=="rs7177338"] <- "Unclassified"


df$Classification <- factor(df$Classification, levels = c("Maternal", "Fetal & Maternal", "Unclassified"))

df1 <- df %>% filter(Classification=="Maternal"|Classification=="Fetal & Maternal")

df2 <- df %>% filter(Classification=="Unclassified")


df1$rsid <- factor(df1$rsid, levels=c("rs72804545","rs2168101","rs180435","rs303998",
                                    "rs11708067","rs138715366","rs2237892"))

df1$gene <- factor(df1$gene, levels=c("EBF1","LMO1","SLC38A4","NLRP13",
                                    "ADCY5","YKT6","KCNQ1"))
  
  

#### Order data for plots

df1 <- df1 %>% arrange(rsid, Genome, Analysis)

df1$Person <- factor(df1$Person, levels = rev(unique(df1$Person)))


#### Italicize the gene names from "gene"


levels(df1$gene) <- paste0("italic('", levels(df1$gene),"')")


##### Plot data


sp1 <- ggplot(df1, aes(Beta, Person)) 

sp1 <- sp1 + geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill=Classification)) 

sp1 <- sp1 + aes(x=Beta, xmin=Beta-1.96*se, xmax=Beta+1.96*se, y=Person, color=Genome) + geom_pointrange(aes(shape=Analysis), size = 1.5) 


sp1 <- sp1 + facet_grid(rsid+gene~.,scales="free",space="free",labeller= label_parsed)

sp1 <- sp1 + geom_vline(xintercept=0) # add 0 line

sp1 <- sp1 + xlab("Effect Size") + ylab("") + scale_color_manual(values=c("#009E73","#D55E00","#0072B2")) + 
  scale_x_continuous(breaks=seq(-0.2,0.4,0.1)) +
  theme_bw() + geom_vline(xintercept = seq(-0.2,0.4,0.1), color ="gray30", linetype = "dashed") +
  theme(text = element_text(face = "bold", size=36, color = "black")) +
  theme(strip.text.y.right = element_text(angle=90, size=28, face="bold")) +
  theme(axis.text.y.left = element_blank()) + 
  theme(axis.ticks.y.left = element_blank()) +
  scale_fill_manual(values = c("Maternal" = "linen","Fetal & Maternal" = "grey90","Unclassified" = "lavenderblush1")) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) +
  guides(color = guide_legend(override.aes=list(fill = NA, linetype = 0, size=1)),
         shape = guide_legend(override.aes = list(linetype = 0)), 
         fill = guide_legend(override.aes=list(linetype = c(1,1,1)))) +
  theme(legend.key.size = unit(1.5, 'cm')) +
  theme(legend.text = element_text(size=24))




# save the plot


png("meta/conditional_analysis/Forest_plots/Plots/maternal_mat-fet.png",
    width=900, height=1600)
sp1
dev.off()



##### Unclassified plot


df2$rsid <- factor(df2$rsid, levels=c("rs150138294","rs10925945","rs4953353","rs74457440",
                                     "rs75512885","rs9800506","rs12543725","rs7177338"))

df2$gene <- factor(df2$gene, levels=c("DCST2","CHRM3","EPAS1","PDLIM5",
                                    "EBF1","FKBP5","SLC45A4","FES"))


#### Order data for plots

df2 <- df2 %>% arrange(rsid, Genome, Analysis)

df2$Person <- factor(df2$Person, levels = rev(unique(df2$Person)))


#### Italicize the gene names from "gene"


levels(df2$gene) <- paste0("italic('", levels(df2$gene),"')")


##### Plot data


sp2 <- ggplot(df2, aes(Beta, Person)) 

sp2 <- sp2 + geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill=Classification)) 

sp2 <- sp2 + aes(x=Beta, xmin=Beta-1.96*se, xmax=Beta+1.96*se, y=Person, color=Genome) + geom_pointrange(aes(shape=Analysis), size = 1.5) 


sp2 <- sp2 + facet_grid(rsid+gene~.,scales="free",space="free",labeller= label_parsed)

sp2 <- sp2 + geom_vline(xintercept=0) # add 0 line

sp2 <- sp2 + xlab("Effect Size") + ylab("") + scale_color_manual(values=c("#009E73","#D55E00","#0072B2")) + 
  scale_x_continuous(breaks=seq(-0.2,0.4,0.1)) +
  theme_bw() + geom_vline(xintercept = seq(-0.2,0.4,0.1), color ="gray30", linetype = "dashed") +
  theme(text = element_text(face = "bold", size=36, color = "black")) +
  theme(strip.text.y.right = element_text(angle=90, size=28, face="bold")) +
  theme(axis.text.y.left = element_blank()) + 
  theme(axis.ticks.y.left = element_blank()) +
  scale_fill_manual(values = c("Maternal" = "linen","Fetal & Maternal" = "grey90","Unclassified" = "lavenderblush1")) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) +
  guides(color = guide_legend(override.aes=list(fill = NA, linetype = 0, size=1)),
         shape = guide_legend(override.aes = list(linetype = 0)), 
         fill = guide_legend(override.aes=list(linetype = c(1,1,1)))) +
  theme(legend.key.size = unit(1.5, 'cm')) +
  theme(legend.text = element_text(size=26))



# save the plot


png("meta/conditional_analysis/Forest_plots/Plots/unclass.png",
    width=900, height=1600)
sp2
dev.off()


sp2 <- sp2 + theme(legend.box = "horizontal")

leg <- as_ggplot(get_legend(sp2))

sp1 <- sp1 + theme(legend.position="none")
sp2 <- sp2 + theme(legend.position="none")



layout <- c(area(1,1),
            area(2,1,8,1),
            area(1,2,8,2)
)

tf <- (leg/sp1) + sp2 + plot_layout(design=layout)

png("Classification_Forest_plots/Plots/SuppFig_3B_Classifications of remaining Loci.png",
    res=300, width=55, height=65, units="cm")
tf
dev.off()

















