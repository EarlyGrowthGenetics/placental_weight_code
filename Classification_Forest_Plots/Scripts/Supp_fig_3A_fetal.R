
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


dff1 <- filter(df, rsid=="rs723177"|rsid=="rs1655296"|rsid=="rs9817452"|rsid=="rs7722058"|rsid=="rs72801474"|rsid=="rs3822394"|
                   rsid=="rs67265526"|rsid=="rs12529634")

dff2 <- filter(df, rsid=="rs11756568"|rsid=="rs1021508"|rsid=="rs6456014"|rsid=="rs7783810"|rsid=="rs10486660"|rsid=="rs6557677"|
                   rsid=="rs1434836"|rsid=="rs1801253")

dff3 <- filter(df, rsid=="rs55958435"|rsid=="rs112635299"|rsid=="rs57790054"|rsid=="rs11866404"|rsid=="rs876987"|
                   rsid=="rs6040436"|rsid=="rs2207099"|rsid=="rs6078190")

nudt <- filter(df, rsid=="rs541641049"|rsid=="rs140691414")

tsnax <- filter(df, rsid=="rs140691414")


dfm <- filter(df, rsid=="rs72804545"|rsid=="rs2168101"|rsid=="rs180435"|rsid=="rs303998")

dfmf <- filter(df, rsid=="rs11708067"|rsid=="rs138715366"|rsid=="rs2237892")


dfu <- filter(df, rsid=="rs150138294"|rsid=="rs10925945"|rsid=="rs4953353"|rsid=="rs74457440"|rsid=="rs75512885"|rsid=="rs9800506"|
                  rsid=="rs12543725"|rsid=="rs7177338")



############## Fetal Plots - first 13

#### Create column for analysis type in plots

dff1$Analysis <- ifelse(grepl('WLM', dff1$Person), "WLM", "Meta")

dff1$Analysis[dff1$Analysis=="Meta"] <- "Meta-analysis"
dff1$Analysis[dff1$Analysis=="WLM"] <- "WLM-adjusted"


#### Create column for facet in plots

dff1$Genome[(grepl('Fetal', dff1$Person))] <- "Fetal"
dff1$Genome[(grepl('Maternal', dff1$Person))] <- "Maternal"
dff1$Genome[(grepl('Paternal', dff1$Person))] <- "Paternal"



dff1$rsid <- factor(dff1$rsid, levels=c("rs723177","rs1655296","rs9817452",
                                        "rs7722058","rs72801474","rs3822394","rs67265526","rs12529634",
                                        "rs11756568","rs1021508","rs6456014","rs7783810"))


dff1$gene[dff1$rsid=="rs723177"] <-"RPL31P11"
dff1$gene[dff1$rsid=="rs1655296"] <- "TSNAX-DISC1"
dff1$gene[dff1$rsid=="rs9817452"] <- "LOC339894"
dff1$gene[dff1$rsid=="rs7722058"] <- "ACTBL2"
dff1$gene[dff1$rsid=="rs72801474"] <- "HSPA4"
dff1$gene[dff1$rsid=="rs3822394"] <- "ARHGAP26"
dff1$gene[dff1$rsid=="rs67265526"] <- "EBF1"
dff1$gene[dff1$rsid=="rs12529634"] <- "HACE1"
dff1$gene[dff1$rsid=="rs11756568"] <- "ESR1"
dff1$gene[dff1$rsid=="rs1021508"] <- "PDE10A"
dff1$gene[dff1$rsid=="rs6456014"] <- "PDE10A"
dff1$gene[dff1$rsid=="rs7783810"] <- "ISPD"


dff1$gene <- factor(dff1$gene, levels=c("RPL31P11","TSNAX-DISC1","LOC339894",
                                        "ACTBL2","HSPA4","ARHGAP26","EBF1","HACE1",
                                        "ESR1","PDE10A","ISPD"))

#### Order data for plots

dff1 <- dff1 %>% arrange(gene, Genome, Analysis)

dff1$Person <- factor(dff1$Person, levels = rev(unique(dff1$Person)))


#### Italicize the gene names from "gene"


levels(dff1$gene) <- paste0("italic('", levels(dff1$gene),"')")

###### Plot data


spf1 <- ggplot(dff1, aes(Beta, Person)) 

spf1 <- spf1 + aes(x=Beta, xmin=Beta-1.96*se, xmax=Beta+1.96*se, y=Person, color=Genome) + geom_pointrange(aes(shape=Analysis), size = 1.5)  +
                xlim(-0.08,0.225) 
            
spf1 <- spf1 + facet_grid(rsid+gene~.,scales="free",space="free",labeller= label_parsed)


spf1 <- spf1 + geom_vline(xintercept=0) # add 0 line

spf1 <- spf1 + xlab("") + ylab("") + 
        theme_bw() + theme(panel.background = element_rect(fill = "honeydew2"), 
                           panel.grid.major.x=element_line(color ="gray30", linetype = "dashed")) +
        scale_color_manual(values=c("#009E73","#D55E00","#0072B2")) + 
        theme(text = element_text(face = "bold", size=32, color = "black")) +
        theme(strip.text.y.right = element_text(angle=90, size=24)) +
        theme(axis.text.y.left = element_blank()) + 
        theme(axis.ticks.y.left = element_blank()) +
        theme(legend.key = element_rect(fill = "white", colour = "black")) +
        guides(color = guide_legend(override.aes=list(fill = NA, linetype = 0, size=1)),
               shape = guide_legend(override.aes = list(linetype = 0)), 
               fill = guide_legend(override.aes=list(linetype = c(1,1,1)))) +
        theme(legend.key.size = unit(1.5, 'cm')) +
        theme(legend.text = element_text(size=20))
  



# save the plot


png("meta/conditional_analysis/Forest_plots/Plots/fetal_1.png",
    width=900, height=1500)
spf1
dev.off()



############## Fetal Plots - second 13

#### Create column for analysis type in plots

dff2$Analysis <- ifelse(grepl('WLM', dff2$Person), "WLM", "Meta")

dff2$Analysis[dff2$Analysis=="Meta"] <- "Meta-analysis"
dff2$Analysis[dff2$Analysis=="WLM"] <- "WLM-adjusted"

#### Create column for facet in plots

dff2$Genome[(grepl('Fetal', dff2$Person))] <- "Fetal"
dff2$Genome[(grepl('Maternal', dff2$Person))] <- "Maternal"
dff2$Genome[(grepl('Paternal', dff2$Person))] <- "Paternal"


dff2$rsid <- factor(dff2$rsid, levels=c("rs11756568","rs1021508","rs6456014","rs7783810",
                                        "rs10486660","rs6557677","rs1434836","rs1801253"))


dff2$gene[dff2$rsid=="rs11756568"] <- "ESR1"
dff2$gene[dff2$rsid=="rs1021508"] <- "PDE10A"
dff2$gene[dff2$rsid=="rs6456014"] <- "PDE10A"
dff2$gene[dff2$rsid=="rs7783810"] <- "ISPD"
dff2$gene[dff2$rsid=="rs10486660"] <- "TBX20"
dff2$gene[dff2$rsid=="rs6557677"] <- "ENTPD4"
dff2$gene[dff2$rsid=="rs1434836"] <- "KLF4"
dff2$gene[dff2$rsid=="rs1801253"] <- "ADRB1"


dff2$gene <- factor(dff2$gene, levels=c("ESR1","PDE10A","ISPD",
                                        "TBX20","ENTPD4","KLF4","ADRB1"))

#### Order data for plots

dff2 <- dff2 %>% arrange(gene, Genome, Analysis)

dff2$Person <- factor(dff2$Person, levels = rev(unique(dff2$Person)))


#### Italicize the gene names from "gene"


levels(dff2$gene) <- paste0("italic('", levels(dff2$gene),"')")

###### Plot data


spf2 <- ggplot(dff2, aes(Beta, Person)) 


spf2 <- spf2 + aes(x=Beta, xmin=Beta-1.96*se, xmax=Beta+1.96*se, y=Person, color=Genome) + geom_pointrange(aes(shape=Analysis), size = 1.5) +
               xlim(-0.08,0.225) 

spf2 <- spf2 + facet_grid(rsid+gene~.,scales="free",space="free",labeller= label_parsed)

spf2 <- spf2 + geom_vline(xintercept=0) # add 0 line

spf2 <- spf2 + xlab("Effect Size") + ylab("") + 
  theme_bw() + theme(panel.background = element_rect(fill = "honeydew2"), 
                     panel.grid.major.x=element_line(color ="gray30", linetype = "dashed")) +
  scale_color_manual(values=c("#009E73","#D55E00","#0072B2")) + 
  theme(text = element_text(face = "bold", size=32, color = "black")) +
  theme(strip.text.y.right = element_text(angle=90, size=24)) +
  theme(axis.text.y.left = element_blank()) + 
  theme(axis.ticks.y.left = element_blank()) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) +
  guides(color = guide_legend(override.aes=list(fill = NA, linetype = 0, size=1)),
         shape = guide_legend(override.aes = list(linetype = 0)), 
         fill = guide_legend(override.aes=list(linetype = c(1,1,1)))) +
  theme(legend.key.size = unit(1.5, 'cm')) +
  theme(legend.text = element_text(size=20))

spf2


# save the plot


png("meta/conditional_analysis/Forest_plots/Plots/fetal_2.png",
    width=900, height=1500)
spf2
dev.off()


############## Fetal Plots - third 8

#### Create column for analysis type in plots

dff3$Analysis <- ifelse(grepl('WLM', dff3$Person), "WLM", "Meta")

dff3$Analysis[dff3$Analysis=="Meta"] <- "Meta-analysis"
dff3$Analysis[dff3$Analysis=="WLM"] <- "WLM-adjusted"

#### Create column for facet in plots

dff3$Genome[(grepl('Fetal', dff3$Person))] <- "Fetal"
dff3$Genome[(grepl('Maternal', dff3$Person))] <- "Maternal"
dff3$Genome[(grepl('Paternal', dff3$Person))] <- "Paternal"



dff3$rsid <- factor(dff3$rsid, levels=c("rs55958435","rs112635299","rs57790054","rs11866404",
                                        "rs876987","rs6040436","rs2207099","rs6078190"))


dff3$gene[dff3$rsid=="rs55958435"] <- "NR2F2"
dff3$gene[dff3$rsid=="rs112635299"] <- "SERPINA1"
dff3$gene[dff3$rsid=="rs57790054"] <- "GPR139"
dff3$gene[dff3$rsid=="rs11866404"] <- "SLC6A2"
dff3$gene[dff3$rsid=="rs876987"] <- "SLC7A5"
dff3$gene[dff3$rsid=="rs6040436"] <- "LOC339593"
dff3$gene[dff3$rsid=="rs2207099"] <- "LOC339593"
dff3$gene[dff3$rsid=="rs6078190"] <- "LOC339593"


dff3$gene <- factor(dff3$gene, levels=c("NR2F2","SERPINA1","GPR139","SLC6A2",
                                        "SLC7A5","LOC339593"))

#### Order data for plots

dff3 <- dff3 %>% arrange(gene, Genome, Analysis)

dff3$Person <- factor(dff3$Person, levels = rev(unique(dff3$Person)))


#### Italicize the gene names from "gene"


levels(dff3$gene) <- paste0("italic('", levels(dff3$gene),"')")

###### Plot data


spf3 <- ggplot(dff3, aes(Beta, Person)) 


spf3 <- spf3 + aes(x=Beta, xmin=Beta-1.96*se, xmax=Beta+1.96*se, y=Person, color=Genome) + geom_pointrange(aes(shape=Analysis), size = 1.5) +
  xlim(-0.08,0.225) 


spf3 <- spf3 + facet_grid(rsid+gene~.,scales="free",space="free",labeller= label_parsed)

spf3 <- spf3 + geom_vline(xintercept=0) # add 0 line

spf3 <- spf3 + xlab("Effect Size") + ylab("") + 
  theme_bw() + theme(panel.background = element_rect(fill = "honeydew2"), 
                     panel.grid.major.x=element_line(color ="gray30", linetype = "dashed")) +
  scale_color_manual(values=c("#009E73","#D55E00","#0072B2")) + 
  theme(text = element_text(face = "bold", size=32, color = "black")) +
  theme(strip.text.y.right = element_text(angle=90, size=24)) +
  theme(axis.text.y.left = element_blank()) + 
  theme(axis.ticks.y.left = element_blank()) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) +
  guides(color = guide_legend(override.aes=list(fill = NA, linetype = 0, size=1)),
         shape = guide_legend(override.aes = list(linetype = 0)), 
         fill = guide_legend(override.aes=list(linetype = c(1,1,1)))) +
  theme(legend.key.size = unit(1.5, 'cm')) +
  theme(legend.text = element_text(size=20))

spf3


# save the plot


png("meta/conditional_analysis/Forest_plots/Plots/fetal_3.png",
    width=900, height=1500)
spf3
dev.off()






#### Create plots for NUDT3 which is on a different scale


#### Create column for analysis type in plots

nudt$Analysis <- ifelse(grepl('WLM', nudt$Person), "WLM", "Meta")

nudt$Analysis[nudt$Analysis=="Meta"] <- "Meta-analysis"
nudt$Analysis[nudt$Analysis=="WLM"] <- "WLM-adjusted"

#### Create column for facet in plots

nudt$Genome[(grepl('Fetal', nudt$Person))] <- "Fetal"
nudt$Genome[(grepl('Maternal', nudt$Person))] <- "Maternal"
nudt$Genome[(grepl('Paternal', nudt$Person))] <- "Paternal"


nudt$rsid <- factor(nudt$rsid, levels=c("rs140691414","rs541641049",))

nudt$gene[nudt$rsid=="rs541641049"] <- "NUDT3"
nudt$gene[nudt$rsid=="rs140691414"] <- "TSNAX-DISC1"

nudt$gene <- factor(nudt$gene, levels=c("TSNAX-DISC1","NUDT3"))

#### Order data for plots

nudt <- nudt %>% arrange(gene, Genome, Analysis)

nudt$Person <- factor(nudt$Person, levels = rev(unique(nudt$Person)))


#### Italicize the gene names from "gene"


levels(nudt$gene) <- paste0("italic('", levels(nudt$gene),"')")

#### Plot data


spf4 <- ggplot(nudt, aes(Beta, Person)) 

spf4 <- spf4 + aes(x=Beta, xmin=Beta-1.96*se, xmax=Beta+1.96*se, y=Person, color=Genome) + geom_pointrange(aes(shape=Analysis), size = 1.5) 


spf4 <- spf4 + facet_grid(rsid+gene~.,scales="free",space="free",labeller= label_parsed)

spf4 <- spf4 + geom_vline(xintercept=0,lty=2) # add 0 line

spf4 <- spf4 + xlab("Effect Size") + ylab("") + 
  theme_bw() + theme(panel.background = element_rect(fill = "honeydew2"),
                     panel.grid.major.x=element_line(color ="gray30", linetype = "dashed")) +
  scale_color_manual(values=c("#009E73","#D55E00","#0072B2")) + 
  theme(text = element_text(face = "bold", size=32, color = "black")) +
  theme(strip.text.y.right = element_text(angle=90, size=24)) +
  theme(axis.text.y.left = element_blank()) + 
  theme(axis.ticks.y.left = element_blank()) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) +
  guides(color = guide_legend(override.aes=list(fill = NA, linetype = 0, size=1)),
         shape = guide_legend(override.aes = list(linetype = 0)), 
         fill = guide_legend(override.aes=list(linetype = c(1,1,1)))) +
  theme(legend.key.size = unit(1.5, 'cm')) +
  theme(legend.text = element_text(size=36))

spf4 <- spf4 + theme(legend.box = "horizontal")

leg <- as_ggplot(get_legend(spf4))

spf4 <- spf4 + theme(legend.position="none")

# save the plot


png("meta/conditional_analysis/Forest_plots/Plots/fetal_nudt3.png",
     width=900, height=300)
spf4
dev.off()




#### Combine the two plots for fetal

spf1 <- spf1 + theme(legend.position="none")

spf2 <- spf2 + theme(legend.position="none")

spf3 <- spf3 + theme(legend.position="none")

tf <- (spf1 + spf2 + spf3) / (spf4 + plot_spacer() + leg) + plot_layout(heights = unit(c(50, 13), c('cm')) )

png("Classification_Forest_plots/Plots/SuppFig_3A_Fetal Classified Loci.png",
    res=300, width=55, height=70, units="cm")
tf
dev.off()
