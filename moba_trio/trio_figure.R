
##
#
# This script builds the figure for the trio analysis.
#
##


# Libraries

library(conflicted)
library(glue)
library(janitor)
library(rjson)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggtext)
library(scico)
library(RColorBrewer)
library(grid)
library(patchwork)


# Namespace conflits

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")


# Load data

trioDF <- read.table(
    file = "docs/figures/trio/merge_decode_06.07.22.txt",
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>% 
    clean_names()

metaAllelesDF <- read.table(
    file = "docs/figures/trio/alleles_meta.txt",
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>% 
    clean_names()


# Remove deCode data for bad proxies

ldThreshold <- 0.2

trioDF$decode_mt_beta[trioDF$proxy_decode_r2 < ldThreshold] <- NA
trioDF$decode_ft_beta[trioDF$proxy_decode_r2 < ldThreshold] <- NA
trioDF$decode_mnt_beta[trioDF$proxy_decode_r2 < ldThreshold] <- NA
trioDF$decode_mt_se[trioDF$proxy_decode_r2 < ldThreshold] <- NA
trioDF$decode_ft_se[trioDF$proxy_decode_r2 < ldThreshold] <- NA
trioDF$decode_mnt_se[trioDF$proxy_decode_r2 < ldThreshold] <- NA
trioDF$decode_mt_p[trioDF$proxy_decode_r2 < ldThreshold] <- NA
trioDF$decode_ft_p[trioDF$proxy_decode_r2 < ldThreshold] <- NA
trioDF$decode_mnt_p[trioDF$proxy_decode_r2 < ldThreshold] <- NA


# Cluster by effect sizes

threshold <- 0.05 / (4 * nrow(trioDF))

trioDF$category <- NA

getCategoryMnT <- function(
    n,
    hMnt,
    hMt,
    hFt,
    hFnt,
    hMntSe,
    hMtSe,
    hFtSe,
    hFntSe,
    hMntP,
    hMtP,
    hFtP,
    hFntP
) {
    
    if (sign(hMnt) == sign(hMt)) {
        
        pTemp <- 2 * pt(q = -abs((hMt - hMnt)/hMtSe), df = n - 3)
        
        if (pTemp < threshold) {
            
            if (hFtP < threshold) {
                
                if (hFntP < threshold) {
                    
                    return("MnT_FT_FnT")
                    
                } else {
                    
                    return("MnT_FT")
                    
                }
                
            } else {
                
                if (hFntP < threshold) {
                    
                    return("MnT_FnT")
                    
                } else {
                    
                    return("MnT")
                    
                }
            }
            
        } else {
            
            if (hFtP < threshold) {
                
                if (hFntP < threshold) {
                    
                    return("MnT_MT_FT_FnT")
                    
                } else {
                    
                    return("MnT_MT_FT")
                    
                }
                
            } else {
                
                if (hFntP < threshold) {
                    
                    return("MnT_MT_FnT")
                    
                } else {
                    
                    return("MnT_MT")
                    
                }
            }
        }
    } else {
        if (hMtP < threshold) {
            
            if (hFtP < threshold) {
                
                if (hFntP < threshold) {
                    
                    return("MnT_MT_FT_FnT")
                    
                } else {
                    
                    return("MnT_MT_FT")
                    
                }
                
            } else {
                
                if (hFntP < threshold) {
                    
                    return("MnT_MT_FnT")
                    
                } else {
                    
                    return("MnT_MT")
                    
                }
            }
            
        } else {
            
            if (hFtP < threshold) {
                
                if (hFntP < threshold) {
                    
                    return("MnT_FT_FnT")
                    
                } else {
                    
                    return("MnT_FT")
                    
                }
                
            } else {
                
                if (hFntP < threshold) {
                    
                    return("MnT_FnT")
                    
                } else {
                    
                    return("MnT")
                    
                }
            }
        }
    }
    
    stop("no category found")
}



getCategoryMT <- function(
    n,
    hMnt,
    hMt,
    hFt,
    hFnt,
    hMntSe,
    hMtSe,
    hFtSe,
    hFntSe,
    hMntP,
    hMtP,
    hFtP,
    hFntP
) {
    
    if ((hFtP < hMntP || sign(hMnt) != sign(hMt)) && sign(hFt) == sign(hMt)) {
        
        pTemp <- 2 * pt(q = -abs((hMt - hFt)/hFtSe), df = n - 3)
        
        if (pTemp < threshold) {
            
            pTemp <- 2 * pt(q = -abs((hMt - hMnt)/hMtSe), df = n - 3)
            
            if (pTemp < threshold) {
                
                if (hFntP < threshold) {
                    
                    return("MT_FnT")
                    
                } else {
                    
                    return("MT")
                    
                }
                
            } else {
                
                if (hFntP < threshold) {
                    
                    return("MnT_MT_FnT")
                    
                } else {
                    
                    return("MnT_MT")
                    
                }
            }
            
        } else {
            
            if (hMntP < threshold) {
                
                if (hFntP < threshold) {
                    
                    return("MnT_MT_FT_FnT")
                    
                } else {
                    
                    return("MnT_MT_FT")
                    
                }
                
            } else {
                
                if (hFntP < threshold) {
                    
                    return("MT_FT_FnT")
                    
                } else {
                    
                    return("MT_FT")
                    
                }
            }
        }
    } else if (sign(hMnt) == sign(hMt)) {
        
        pTemp <- 2 * pt(q = -abs((hMt - hMnt)/hMtSe), df = n - 3)
        
        if (pTemp < threshold) {
            
            pTemp <- 2 * pt(q = -abs((hMt - hFt)/hFtSe), df = n - 3)
            
            if (pTemp < threshold) {
                
                if (hFntP < threshold) {
                    
                    return("MT_FnT")
                    
                } else {
                    
                    return("MT")
                    
                }
                
            } else {
                
                if (hFntP < threshold) {
                    
                    return("MT_FT_FnT")
                    
                } else {
                    
                    return("MT_FT")
                    
                }
            }
            
        } else {
            
            if (hFtP < threshold) {
                
                if (hFntP < threshold) {
                    
                    return("MnT_MT_FT_FnT")
                    
                } else {
                    
                    return("MnT_MT_FT")
                    
                }
                
            } else {
                
                if (hFntP < threshold) {
                    
                    return("MnT_MT_FnT")
                    
                } else {
                    
                    return("MnT_MT")
                    
                }
            }
        }
    } else {
        
        if (hFtP < threshold) {
            
            if (hMntP < threshold) {
                
                if (hFntP < threshold) {
                    
                    return("MnT_MT_FT_FnT")
                    
                } else {
                    
                    return("MnT_MT_FT")
                    
                }
                
            } else {
                
                if (hFntP < threshold) {
                    
                    return("MT_FT_FnT")
                    
                } else {
                    
                    return("MT_FT")
                    
                }
            }
        } else {
            
            if (hMntP < threshold) {
                
                if (hFntP < threshold) {
                    
                    return("MnT_MT_FnT")
                    
                } else {
                    
                    return("MnT_MT")
                    
                }
                
            } else {
                
                if (hFntP < threshold) {
                    
                    return("MT_FnT")
                    
                } else {
                    
                    return("MT")
                    
                }
            }
        }
    }
    
    
    stop("no category found")
}



getCategoryFT <- function(
    n,
    hMnt,
    hMt,
    hFt,
    hFnt,
    hMntSe,
    hMtSe,
    hFtSe,
    hFntSe,
    hMntP,
    hMtP,
    hFtP,
    hFntP
) {
    
    if (sign(hFt) == sign(hMt)) {
        
        pTemp <- 2 * pt(q = -abs((hMt - hFt)/hMtSe), df = n - 3)
        
        if (pTemp < threshold) {
            
            if (hMntP < threshold) {
                
                if (hFntP < threshold) {
                    
                    return("MnT_FT_FnT")
                    
                } else {
                    
                    return("MnT_FT")
                    
                }
                
            } else {
                
                if (hFntP < threshold) {
                    
                    return("FT_FnT")
                    
                } else {
                    
                    return("FT")
                    
                }
            }
            
        } else {
            
            if (hMntP < threshold) {
                
                if (hFntP < threshold) {
                    
                    return("MnT_MT_FT_FnT")
                    
                } else {
                    
                    return("MnT_MT_FT")
                    
                }
                
            } else {
                
                if (hFntP < threshold) {
                    
                    return("MT_FT_FnT")
                    
                } else {
                    
                    return("MT_FT")
                    
                }
            }
        }
    } else {
        
        if (hMtP < threshold) {
            
            if (hMntP < threshold) {
                
                if (hFntP < threshold) {
                    
                    return("MnT_MT_FT_FnT")
                    
                } else {
                    
                    return("MnT_MT_FT")
                    
                }
                
            } else {
                
                if (hFntP < threshold) {
                    
                    return("MT_FT_FnT")
                    
                } else {
                    
                    return("MT_FT")
                    
                }
            }
            
        } else {
            
            if (hMntP < threshold) {
                
                if (hFntP < threshold) {
                    
                    return("MnT_FT_FnT")
                    
                } else {
                    
                    return("MnT_FT")
                    
                }
                
            } else {
                
                if (hFntP < threshold) {
                    
                    return("FT_FnT")
                    
                } else {
                    
                    return("FT")
                    
                }
            }
        }
    }
    
    
    stop("no category found")
}



getCategoryFnT <- function(
    n,
    hMnt,
    hMt,
    hFt,
    hFnt,
    hMntSe,
    hMtSe,
    hFtSe,
    hFntSe,
    hMntP,
    hMtP,
    hFtP,
    hFntP
) {
    
    if (sign(hFnt) == sign(hFt)) {
        
        pTemp <- 2 * pt(q = -abs((hFnt - hFt)/hFtSe), df = n - 3)
        
        if (pTemp < threshold) {
            
            if (hMntP < threshold) {
                
                if (hMtP < threshold) {
                    
                    return("MnT_MT_FnT")
                    
                } else {
                    
                    
                    return("MnT_FnT")
                    
                }
                
            } else {
                
                if (hMtP < threshold) {
                    
                    return("MT_FnT")
                    
                } else {
                    
                    return("FnT")
                    
                }
            }
            
        } else {
            
            if (hMntP < threshold) {
                
                if (hMtP < threshold) {
                    
                    return("MnT_MT_FT_FnT")
                    
                } else {
                    
                    
                    return("MT_FT_FnT")
                    
                }
                
            } else {
                
                if (hMtP < threshold) {
                    
                    return("MT_FT_FnT")
                    
                } else {
                    
                    
                    return("MT_FT_FnT")
                    
                }
            }
        }
    } else {
        
        if (hFtP < threshold) {
            
            if (hMntP < threshold) {
                
                if (hMtP < threshold) {
                    
                    return("MnT_MT_FT_FnT")
                    
                } else {
                    
                    
                    return("MnT_FT_FnT")
                    
                }
                
            } else {
                
                if (hMtP < threshold) {
                    
                    return("MT_FT_FnT")
                    
                } else {
                    
                    return("FT_FnT")
                    
                }
            }
            
        } else {
            
            if (hMntP < threshold) {
                
                if (hMtP < threshold) {
                    
                    return("MnT_MT_FnT")
                    
                } else {
                    
                    return("MnT_FnT")
                    
                }
                
            } else {
                
                if (hMtP < threshold) {
                    
                    return("MT_FnT")
                    
                } else {
                    
                    return("FnT")
                    
                }
            }
        }
    }
    
    
    stop("no category found")
}




for (i in 1:nrow(trioDF)) {
    
    snp <- trioDF$rsid[i]
    
    n <- trioDF$n[i]
    
    hMnt <- trioDF$h_bmnt_moba_pw[i]
    hMt <- trioDF$h_bmt_moba_pw[i]
    hFt <- trioDF$h_bft_moba_pw[i]
    hFnt <- trioDF$h_bfnt_moba_pw[i]
    
    hMntSe <- trioDF$h_bmnt_se_moba_pw[i]
    hMtSe <- trioDF$h_bmt_se_moba_pw[i]
    hFtSe <- trioDF$h_bft_se_moba_pw[i]
    hFntSe <- trioDF$h_bfnt_se_moba_pw[i]
    
    hMntP <- trioDF$h_bmnt_p_moba_pw[i]
    hMtP <- trioDF$h_bmt_p_moba_pw[i]
    hFtP <- trioDF$h_bft_p_moba_pw[i]
    hFntP <- trioDF$h_bfnt_p_moba_pw[i]
    
    if (!is.na(hMntP)) {
        
        if (hMntP <= hMtP && hMntP <= hFtP && hMntP <= hFntP) {
            
            trioDF$category[i] <- getCategoryMnT(
                n,
                hMnt,
                hMt,
                hFt,
                hFnt,
                hMntSe,
                hMtSe,
                hFtSe,
                hFntSe,
                hMntP,
                hMtP,
                hFtP,
                hFntP
            )
            
        } else if (hMtP <= hMntP && hMtP <= hFtP && hMtP <= hFntP) {
            
            trioDF$category[i] <- getCategoryMT(
                n,
                hMnt,
                hMt,
                hFt,
                hFnt,
                hMntSe,
                hMtSe,
                hFtSe,
                hFntSe,
                hMntP,
                hMtP,
                hFtP,
                hFntP
            )
            
        } else if (hFtP <= hMntP && hFtP <= hMtP && hFtP <= hFntP) {
            
            trioDF$category[i] <- getCategoryFT(
                n,
                hMnt,
                hMt,
                hFt,
                hFnt,
                hMntSe,
                hMtSe,
                hFtSe,
                hFntSe,
                hMntP,
                hMtP,
                hFtP,
                hFntP
            )
            
        } else if (hFntP <= hMntP && hFntP <= hMtP && hFntP <= hFtP) {
            
            trioDF$category[i] <- getCategoryFnT(
                n,
                hMnt,
                hMt,
                hFt,
                hFnt,
                hMntSe,
                hMtSe,
                hFtSe,
                hFntSe,
                hMntP,
                hMtP,
                hFtP,
                hFntP
            )
        }
    }
}


# Swap alleles trio

betaColumns <- c("h_bmnt_moba_pw", "h_bmt_moba_pw", "h_bft_moba_pw", "h_bfnt_moba_pw", "cmf_bc_moba_pw", "cmf_bm_moba_pw", "cmf_bf_moba_pw", "cmf_mt_bc_moba_pw", "cmf_mt_bm_moba_pw", "cmf_mt_bf_moba_pw", "cmf_mt_bmt_moba_pw",
                 "h_bmnt_moba_bw", "h_bmt_moba_bw", "h_bft_moba_bw", "h_bfnt_moba_bw", "cmf_bc_moba_bw", "cmf_bm_moba_bw", "cmf_bf_moba_bw", "cmf_mt_bc_moba_bw", "cmf_mt_bm_moba_bw", "cmf_mt_bf_moba_bw", "cmf_mt_bmt_moba_bw",
                 "decode_ft_beta", "decode_mt_beta", "decode_mnt_beta")

for (i in 1:nrow(trioDF)) {
    
    rsid <- trioDF$rsid[i]
    
    ea <- metaAllelesDF$ea[metaAllelesDF$rsid == rsid]
    
    if (length(ea) != 1) {
        
        stop(glue("Effect allele for {rsid} not found."))
        
    }
    
    oa <- metaAllelesDF$oa[metaAllelesDF$rsid == rsid]
    
    if (length(oa) != 1) {
        
        stop(glue("Other allele for {rsid} not found."))
        
    }
    
    if (trioDF$tested_allele[i] == oa && trioDF$other_allele[i] == ea) {
        
        swap <- T
        
    } else if (trioDF$tested_allele[i] == ea && trioDF$other_allele[i] == oa) {
        
        swap <- F
        
    } else {
        
        stop(glue("Allele mismatch."))
        
    }
    
    if (swap) {
        
        testedAllele <- trioDF$tested_allele[i]
        trioDF$tested_allele[i] <- trioDF$other_allele[i]
        trioDF$other_allele[i] <- testedAllele
        
        for (betaColumn in betaColumns) {
            
            trioDF[[betaColumn]][i] <- -trioDF[[betaColumn]][i]
            
        }
    }
}


# Sort

trioDF$categoryFactor <- factor(trioDF$category, levels = c("MnT_MT", "MT", "MT_FT", "FT", "FnT"))

trioDF$categoryNameFactor <- trioDF$categoryFactor
levels(trioDF$categoryNameFactor) <- c("Maternal", "Maternal Transmitted", "Fetal", "Paternal Transmitted", "Paternal non-Transmitted")
trioDF$categoryName <- as.character(trioDF$categoryNameFactor)

if (sum(is.na(trioDF$categoryFactor)) > 0) {
    
    stop("Missing or unsupported category")
    
}

trioDF <- trioDF %>% 
    arrange(
        categoryFactor
    ) %>% 
    mutate(
        rank = NA
    )

rank <- 1

for (categoryValue in levels(trioDF$categoryFactor)) {
    
    tempDF <- trioDF %>% 
        filter(
            category == categoryValue
        )
    
    if (categoryValue == "MnT_MT") {
        
        tempDF %>% 
            arrange(
                desc(cmf_bm_moba_pw), cmf_bc_moba_pw
            )
        
    } else if (categoryValue == "MT") {
        
        tempDF %>% 
            arrange(
                desc(h_bmt_moba_pw), cmf_bc_moba_pw
            )
        
    } else if (categoryValue == "MT_FT") {
        
        tempDF %>% 
            arrange(
                cmf_bc_moba_pw
            )
        
    } else if (categoryValue == "FT") {
        
        tempDF %>% 
            arrange(
                h_bft_moba_pw
            )
        
    } else if (categoryValue == "FnT") {
        
        tempDF %>% 
            arrange(
                h_bfnt_moba_pw
            )
        
    } else {
        
        stop("Category not supported")
        
    }
    
    for (i in 1:nrow(tempDF)) {
        
        rsid <- tempDF$rsid[i]
        
        trioDF$rank[which(trioDF$rsid == rsid)] <- rank
        
        rank <- rank + 1
        
    }
}


# Forest of all estimates

trioDF <- trioDF %>% 
    arrange(
        contig, position
    )

dx <- 0.2
margin <- 0.75

plotPatchwork <- NULL
plotPatchwork1 <- NULL
plotPatchwork2 <- NULL
plotPatchwork3 <- NULL

for (i in 1:nrow(trioDF)) {
    
    locus <- trioDF$locus[i]
    rsid <- trioDF$rsid[i]
    maternal_fetal_category <- trioDF$maternal_fetal_classification[i]
    pw_bw_category <- trioDF$pw_bw_classification[i]
    
    # Plot haplotypes
    
    if (!is.na(trioDF$decode_mt_beta[i])) {
        
        if (-min(trioDF$decode_mnt_beta[i], trioDF$decode_mt_beta[i], trioDF$decode_ft_beta[i]) > max(trioDF$decode_mnt_beta[i], trioDF$decode_mt_beta[i], trioDF$decode_ft_beta[i])) {
            
            trioDF$decode_mnt_beta[i] <- -trioDF$decode_mnt_beta[i]
            trioDF$decode_mt_beta[i] <- -trioDF$decode_mt_beta[i]
            trioDF$decode_ft_beta[i] <- -trioDF$decode_ft_beta[i]
            
        }
        
        plotDF <- data.frame(
            x = c(1:4 - dx, 1:4, 1:3 + dx),
            beta = c(
                trioDF$h_bmnt_moba_pw[i],
                trioDF$h_bmt_moba_pw[i],
                trioDF$h_bft_moba_pw[i],
                trioDF$h_bfnt_moba_pw[i],
                trioDF$h_bmnt_moba_bw[i],
                trioDF$h_bmt_moba_bw[i],
                trioDF$h_bft_moba_bw[i],
                trioDF$h_bfnt_moba_bw[i],
                trioDF$decode_mnt_beta[i],
                trioDF$decode_mt_beta[i],
                trioDF$decode_ft_beta[i]
            ),
            se = c(
                trioDF$h_bmnt_se_moba_pw[i],
                trioDF$h_bmt_se_moba_pw[i],
                trioDF$h_bft_se_moba_pw[i],
                trioDF$h_bfnt_se_moba_pw[i],
                trioDF$h_bmnt_se_moba_bw[i],
                trioDF$h_bmt_se_moba_bw[i],
                trioDF$h_bft_se_moba_bw[i],
                trioDF$h_bfnt_se_moba_bw[i],
                trioDF$decode_mnt_se[i],
                trioDF$decode_mt_se[i],
                trioDF$decode_ft_se[i]
            ),
            col = c(
                rep("Placenta Weight (MoBa)", 4),
                rep("Birth Weight (MoBa)", 4),
                rep("Birth Weight (deCode)", 3)
            ),
            stringsAsFactors = F
        ) %>% 
            mutate(
                ymin = beta - qnorm(0.975) * se,
                ymax = beta + qnorm(0.975) * se
            )
        
        plotDF$colorFactor <- factor(plotDF$col, levels = c("Placenta Weight (MoBa)", "Birth Weight (MoBa)", "Birth Weight (deCode)"))
        
        colors <- c("darkblue", "grey40", "grey60")
        
    } else {
        
        plotDF <- data.frame(
            x = c(1:4 - dx/2, 1:4 + dx/2),
            beta = c(
                trioDF$h_bmnt_moba_pw[i],
                trioDF$h_bmt_moba_pw[i],
                trioDF$h_bft_moba_pw[i],
                trioDF$h_bfnt_moba_pw[i],
                trioDF$h_bmnt_moba_bw[i],
                trioDF$h_bmt_moba_bw[i],
                trioDF$h_bft_moba_bw[i],
                trioDF$h_bfnt_moba_bw[i]
            ),
            se = c(
                trioDF$h_bmnt_se_moba_pw[i],
                trioDF$h_bmt_se_moba_pw[i],
                trioDF$h_bft_se_moba_pw[i],
                trioDF$h_bfnt_se_moba_pw[i],
                trioDF$h_bmnt_se_moba_bw[i],
                trioDF$h_bmt_se_moba_bw[i],
                trioDF$h_bft_se_moba_bw[i],
                trioDF$h_bfnt_se_moba_bw[i]
            ),
            col = c(
                rep("Placenta Weight (MoBa)", 4),
                rep("Birth Weight (MoBa)", 4)
            ),
            stringsAsFactors = F
        ) %>% 
            mutate(
                ymin = beta - qnorm(0.975) * se,
                ymax = beta + qnorm(0.975) * se
            )
        
        plotDF$colorFactor <- factor(plotDF$col, levels = c("Placenta Weight (MoBa)", "Birth Weight (MoBa)"))
        
        colors <- c("darkblue", "grey40")
        
    }
    
    plot <- ggplot() +
        theme_bw() +
        geom_hline(
            yintercept = 0
        ) +
        geom_segment(
            data = plotDF,
            mapping = aes(
                x = x,
                xend = x,
                y = ymin,
                yend = ymax,
                col = colorFactor
            )
        ) +
        geom_point(
            data = plotDF,
            mapping = aes(
                x = x,
                y = beta,
                col = colorFactor
            )
        ) +
        scale_x_continuous(
            breaks = 1:4,
            labels = c("Maternal\nnon-Transmitted", "Maternal\nTransmitted", "Paternal\nTransmitted", "Paternal\nnon-Transmitted")
        ) +
        scale_y_continuous(
            name = "Beta [95% CI]"
        ) +
        scale_color_manual(
            values = colors
        ) + 
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            axis.title.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.border = element_blank(),
            plot.margin = margin(t = margin, r = margin, b = 1.5 * margin, l = margin, unit = "cm")
        ) + 
        labs(
            title = glue("{locus} ({rsid}) - Haplotypes"),
            subtitle = glue("{maternal_fetal_category}, {pw_bw_category}")
        )
    
    if (!is.null(plotPatchwork)) {
        
        plotPatchwork <- plotPatchwork + plot 
        
    } else {
        
        plotPatchwork <- plot
        
    }
    
    if (i <= 16) {
        
        if (!is.null(plotPatchwork1)) {
            
            plotPatchwork1 <- plotPatchwork1 + plot 
            
        } else {
            
            plotPatchwork1 <- plot
            
        }
        
    } else if (i <= 32) {
        
        if (!is.null(plotPatchwork2)) {
            
            plotPatchwork2 <- plotPatchwork2 + plot 
            
        } else {
            
            plotPatchwork2 <- plot
            
        }
        
    } else {
        
        if (!is.null(plotPatchwork3)) {
            
            plotPatchwork3 <- plotPatchwork3 + plot 
            
        } else {
            
            plotPatchwork3 <- plot
            
        }
        
    }
    
    
    # Plot trio + mt
    
    plotDF <- data.frame(
        x = c(1:4 - dx/2, 1:4 + dx/2),
        beta = c(
            trioDF$cmf_mt_bm_moba_pw[i],
            trioDF$cmf_mt_bc_moba_pw[i],
            trioDF$cmf_mt_bf_moba_pw[i],
            trioDF$cmf_mt_bmt_moba_pw[i],
            trioDF$cmf_mt_bm_moba_bw[i],
            trioDF$cmf_mt_bc_moba_bw[i],
            trioDF$cmf_mt_bf_moba_bw[i],
            trioDF$cmf_mt_bmt_moba_bw[i]
        ),
        se = c(
            trioDF$cmf_mt_bm_se_moba_pw[i],
            trioDF$cmf_mt_bc_se_moba_pw[i],
            trioDF$cmf_mt_bf_se_moba_pw[i],
            trioDF$cmf_mt_bmt_se_moba_pw[i],
            trioDF$cmf_mt_bm_se_moba_bw[i],
            trioDF$cmf_mt_bc_se_moba_bw[i],
            trioDF$cmf_mt_bf_se_moba_bw[i],
            trioDF$cmf_mt_bmt_se_moba_bw[i]
        ),
        col = c(
            rep("Placenta Weight (MoBa)", 4),
            rep("Birth Weight (MoBa)", 4)
        ),
        stringsAsFactors = F
    ) %>% 
        mutate(
            ymin = beta - qnorm(0.975) * se,
            ymax = beta + qnorm(0.975) * se
        )
    
    plotDF$colorFactor <- factor(plotDF$col, levels = c("Placenta Weight (MoBa)", "Birth Weight (MoBa)"))
    
    colors <- c("darkblue", "grey40")
    
    plot <- ggplot() +
        theme_bw() +
        geom_hline(
            yintercept = 0
        ) +
        geom_segment(
            data = plotDF,
            mapping = aes(
                x = x,
                xend = x,
                y = ymin,
                yend = ymax,
                col = colorFactor
            )
        ) +
        geom_point(
            data = plotDF,
            mapping = aes(
                x = x,
                y = beta,
                col = colorFactor
            )
        ) +
        scale_x_continuous(
            breaks = 1:4,
            labels = c("Mother", "Child", "Father", "Maternal\nTransmitted")
        ) +
        scale_y_continuous(
            name = "Beta [95% CI]"
        ) +
        scale_color_manual(
            values = colors
        ) + 
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            axis.title.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.border = element_blank(),
            plot.margin = margin(t = margin, r = margin, b = 1.5 * margin, l = 0, unit = "cm")
        ) + 
        labs(
            title = glue("Trio + imprinting")
        )
    
    plotPatchwork <- plotPatchwork + plot
    
    if (i <= 16) {
        
        plotPatchwork1 <- plotPatchwork1 + plot
        
    } else if (i <= 32) {
        
        plotPatchwork2 <- plotPatchwork2 + plot
        
    } else {
        
        plotPatchwork3 <- plotPatchwork3 + plot
        
    }
    
}

plotPatchwork <- plotPatchwork + 
    plot_layout(
        ncol = 4
    )
plotPatchwork1 <- plotPatchwork1 + 
    plot_layout(
        ncol = 4
    )
plotPatchwork2 <- plotPatchwork2 + 
    plot_layout(
        ncol = 4
    )
plotPatchwork3 <- plotPatchwork3 + 
    plot_layout(
        ncol = 4
    )


pdf(
    file = glue("docs/figures/trio/haplotypes_supplementary_03.02.22.pdf"),
    width = unit(21, "cm"),
    height = unit(3 * 29.7, "cm"), 
    pointsize = 10
)
grid.draw(plotPatchwork)
dummy <- dev.off()

pdf(
    file = glue("docs/figures/trio/haplotypes_supplementary_1_03.02.22.pdf"),
    width = unit(21, "cm"),
    height = unit(29.7, "cm"), 
    pointsize = 10
)
grid.draw(plotPatchwork1)
dummy <- dev.off()
pdf(
    file = glue("docs/figures/trio/haplotypes_supplementary_2_03.02.22.pdf"),
    width = unit(21, "cm"),
    height = unit(29.7, "cm"), 
    pointsize = 10
)
grid.draw(plotPatchwork2)
dummy <- dev.off()
pdf(
    file = glue("docs/figures/trio/haplotypes_supplementary_3_03.02.22.pdf"),
    width = unit(21, "cm"),
    height = unit(29.7*5/8, "cm"), 
    pointsize = 10
)
grid.draw(plotPatchwork3)
dummy <- dev.off()



# Figure for the main text

trioDF <- trioDF %>% 
    mutate(
        maternal_fetal_classification_factor = factor(maternal_fetal_classification, levels = c("Maternal", "Fetal & Maternal", "Fetal & Maternal opposite directions", "Fetal", "Unclassified")),
        pw_bw_classification = factor(pw_bw_classification, levels = c("Predominantly or only PW", "PW & BW same direction", "PW & BW opposite directions"))
    ) %>% 
    arrange(
        pw_bw_classification, maternal_fetal_classification_factor, rank
    ) %>% 
    mutate(
        rank = row_number() + 0.5 * (as.numeric(pw_bw_classification) - 1)
    ) %>% 
    separate(
        locus, 
        into = c("locus"), 
        extra = "drop"
    )

pw_fm_min <- min(trioDF$rank[(trioDF$maternal_fetal_classification == "Fetal & Maternal" | trioDF$maternal_fetal_classification == "Fetal & Maternal opposite directions") & trioDF$pw_bw_classification == "Predominantly or only PW"])
pw_fm_max <- max(trioDF$rank[(trioDF$maternal_fetal_classification == "Fetal & Maternal" | trioDF$maternal_fetal_classification == "Fetal & Maternal opposite directions") & trioDF$pw_bw_classification == "Predominantly or only PW"])
pw_f_min <- min(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal" & trioDF$pw_bw_classification == "Predominantly or only PW"])
pw_f_max <- max(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal" & trioDF$pw_bw_classification == "Predominantly or only PW"])
pw_u_min <- min(trioDF$rank[trioDF$maternal_fetal_classification == "Unclassified" & trioDF$pw_bw_classification == "Predominantly or only PW"])
pw_u_max <- max(trioDF$rank[trioDF$maternal_fetal_classification == "Unclassified" & trioDF$pw_bw_classification == "Predominantly or only PW"])

pw_bw_m_min <- min(trioDF$rank[trioDF$maternal_fetal_classification == "Maternal" & trioDF$pw_bw_classification == "PW & BW same direction"])
pw_bw_m_max <- max(trioDF$rank[trioDF$maternal_fetal_classification == "Maternal" & trioDF$pw_bw_classification == "PW & BW same direction"])
pw_bw_fm_min <- min(trioDF$rank[(trioDF$maternal_fetal_classification == "Fetal & Maternal" | trioDF$maternal_fetal_classification == "Fetal & Maternal opposite directions") & trioDF$pw_bw_classification == "PW & BW same direction"])
pw_bw_fm_max <- max(trioDF$rank[(trioDF$maternal_fetal_classification == "Fetal & Maternal" | trioDF$maternal_fetal_classification == "Fetal & Maternal opposite directions") & trioDF$pw_bw_classification == "PW & BW same direction"])
pw_bw_f_min <- min(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal" & trioDF$pw_bw_classification == "PW & BW same direction"])
pw_bw_f_max <- max(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal" & trioDF$pw_bw_classification == "PW & BW same direction"])
pw_bw_u_min <- min(trioDF$rank[trioDF$maternal_fetal_classification == "Unclassified" & trioDF$pw_bw_classification == "PW & BW same direction"])
pw_bw_u_max <- max(trioDF$rank[trioDF$maternal_fetal_classification == "Unclassified" & trioDF$pw_bw_classification == "PW & BW same direction"])

pw_bw_o_f_min <- min(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal" & trioDF$pw_bw_classification == "PW & BW opposite directions"])
pw_bw_o_f_max <- max(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal" & trioDF$pw_bw_classification == "PW & BW opposite directions"])


maxRank <- max(trioDF$rank)

betaDF <- data.frame(
    rank = c(trioDF$rank, trioDF$rank, trioDF$rank, trioDF$rank, 
             trioDF$rank, trioDF$rank, trioDF$rank, trioDF$rank,
             trioDF$rank, trioDF$rank, trioDF$rank),
    variable = c(rep("mnt_moba_pw", nrow(trioDF)), rep("mt_moba_pw", nrow(trioDF)), rep("ft_moba_pw", nrow(trioDF)), rep("fnt_moba_pw", nrow(trioDF)), 
                 rep("mnt_moba_bw", nrow(trioDF)), rep("mt_moba_bw", nrow(trioDF)), rep("ft_moba_bw", nrow(trioDF)), rep("fnt_moba_bw", nrow(trioDF)),
                 rep("mnt_decode_bw", nrow(trioDF)), rep("mt_decode_bw", nrow(trioDF)), rep("ft_decode_bw", nrow(trioDF))),
    beta = c(trioDF$h_bmnt_moba_pw, trioDF$h_bmt_moba_pw, trioDF$h_bft_moba_pw, trioDF$h_bfnt_moba_pw, 
             trioDF$h_bmnt_moba_bw, trioDF$h_bmt_moba_bw, trioDF$h_bft_moba_bw, trioDF$h_bfnt_moba_bw,
             trioDF$decode_mnt_beta, trioDF$decode_mt_beta, trioDF$decode_ft_beta),
    se = c(trioDF$h_bmnt_se_moba_pw, trioDF$h_bmt_se_moba_pw, trioDF$h_bft_se_moba_pw, trioDF$h_bfnt_se_moba_pw, 
           trioDF$h_bmnt_se_moba_bw, trioDF$h_bmt_se_moba_bw, trioDF$h_bft_se_moba_bw, trioDF$h_bfnt_se_moba_bw,
           trioDF$decode_mnt_se, trioDF$decode_mt_se, trioDF$decode_ft_se),
    stringsAsFactors = F
) %>% 
    mutate(
        variableFactor = factor(variable, levels = c("mnt_moba_pw", "mt_moba_pw", "ft_moba_pw", "fnt_moba_pw", "mnt_moba_bw", "mt_moba_bw", "ft_moba_bw", "fnt_moba_bw", "mnt_decode_bw", "mt_decode_bw", "ft_decode_bw")),
        x = as.numeric(variableFactor),
        y = -rank,
        t = beta/se
    ) %>% 
    filter(
        !is.na(beta)
    )

maxBetaAbs <- max(abs(betaDF$beta[!is.na(betaDF$beta)]))
maxTAbs <- max(abs(betaDF$t[!is.na(betaDF$t)]))

heatMap <- ggplot() +
    theme_void() +
    geom_text(
        mapping = aes(
            x = 1:11,
            y = 0,
            label = c("MnT", "MT", "PT", "PnT", "MnT", "MT", "PT", "PnT", "MnT", "MT", "PT")
        ),
        size = 3,
        col = "grey40",
        hjust = 0.5,
        vjust = 0
    ) +
    geom_text(
        mapping = aes(
            x = 11.75,
            y = 0,
            label = "WLM"
        ),
        size = 3,
        col = "grey40",
        hjust = 0,
        vjust = 0
    ) +
    geom_richtext(
        data = trioDF,
        mapping = aes(
            x = 0.4,
            y = -rank,
            label = glue("<i>{locus}</i> ({rsid})")
        ),
        size = 3,
        col = "grey30",
        hjust = 1,
        vjust = 0.5,
        fill = NA, 
        label.color = NA
    ) +
    geom_tile(
        data = betaDF,
        mapping = aes(
            x = x,
            y = y,
            fill = t
        )
    ) +
    geom_segment(
        mapping = aes(
            x = c(0.5, 4.5, 8.5, 11.5),
            xend = c(0.5, 4.5, 8.5, 11.5),
            y = -0.25,
            yend = -pw_u_max - 0.5
        ),
        col = "grey70"
    ) +
    geom_segment(
        mapping = aes(
            x = c(0.5, 4.5, 8.5, 11.5),
            xend = c(0.5, 4.5, 8.5, 11.5),
            y = -pw_bw_m_min + 0.5,
            yend = -pw_bw_u_max - 0.5
        ),
        col = "grey70"
    ) +
    geom_segment(
        mapping = aes(
            x = c(0.5, 4.5, 8.5, 11.5),
            xend = c(0.5, 4.5, 8.5, 11.5),
            y = -pw_bw_o_f_min + 0.5,
            yend = -pw_bw_o_f_max - 0.5
        ),
        col = "grey70"
    ) +
    geom_segment(
        mapping = aes(
            x = 0.5,
            xend = 11.6,
            y = c(-pw_fm_min + 0.5, -pw_fm_max - 0.5, -pw_f_min + 0.5, -pw_f_max - 0.5, -pw_u_min + 0.5, -pw_u_max - 0.5, -pw_bw_m_min + 0.5, -pw_bw_m_max - 0.5, -pw_bw_fm_min + 0.5, -pw_bw_fm_max - 0.5, -pw_bw_f_min + 0.5, -pw_bw_f_max - 0.5, -pw_bw_u_min + 0.5, -pw_bw_u_max - 0.5, -pw_bw_o_f_min + 0.5, -pw_bw_o_f_max - 0.5),
            yend = c(-pw_fm_min + 0.5, -pw_fm_max - 0.5, -pw_f_min + 0.5, -pw_f_max - 0.5, -pw_u_min + 0.5, -pw_u_max - 0.5, -pw_bw_m_min + 0.5, -pw_bw_m_max - 0.5, -pw_bw_fm_min + 0.5, -pw_bw_fm_max - 0.5, -pw_bw_f_min + 0.5, -pw_bw_f_max - 0.5, -pw_bw_u_min + 0.5, -pw_bw_u_max - 0.5, -pw_bw_o_f_min + 0.5, -pw_bw_o_f_max - 0.5)
        ),
        col = "grey70"
    ) +
    geom_text(
        mapping = aes(
            x = 11.75,
            y = c(-(pw_fm_min + pw_fm_max) / 2, -(pw_f_min + pw_f_max) / 2, -(pw_u_min + pw_u_max) / 2, -(pw_bw_m_min + pw_bw_m_max) / 2, -(pw_bw_fm_min + pw_bw_fm_max) / 2, -(pw_bw_f_min + pw_bw_f_max) / 2, -(pw_bw_u_min + pw_bw_u_max) / 2, -(pw_bw_o_f_min + pw_bw_o_f_max) / 2),
            label = c("Fetal & Maternal", "Fetal", "Ambiguous", "Maternal", "Fetal & Maternal\nopposite directions", "Fetal", "Ambiguous", "Fetal")
        ),
        size = 3.5,
        col = "grey30",
        hjust = 0,
        vjust = 0.5
    ) +
    geom_text(
        mapping = aes(
            x = -4.1,
            y = c(-(pw_fm_min + pw_u_max) / 2, -(pw_bw_m_min + pw_bw_u_max) / 2, -(pw_bw_o_f_min + pw_bw_o_f_max) / 2),
            label = c("Predominantly or only PW", "PW & BW same direction", "*")
        ),
        size = 4.5,
        col = "grey30",
        hjust = 0.5,
        vjust = 1,
        angle = 90
    ) +
    geom_text(
        mapping = aes(
            x = -4,
            y = -maxRank - 1.5,
            label = "* PW & BW opposite direction"
        ),
        size = 4.1,
        col = "grey30",
        vjust = 1,
        hjust = 0
    ) +
    geom_richtext(
        mapping = aes(
            x = c(2.5, 6.5, 10),
            y = 2,
            label = c("Placenta Weight<br>MoBa", "Birth Weight<br>MoBa", "Birth Weight<br>Juliusdottir <i>et al.</i>")
        ),
        size = 3.5,
        col = "grey30",
        hjust = 0.5,
        vjust = 0,
        fill = NA, 
        label.color = NA
    ) +
    scale_fill_scico(
        name = "Beta / se",
        palette = "vik",
        limits = c(-maxTAbs, maxTAbs)
    ) +
    coord_cartesian(
        xlim = c(-4, 14.7),
        ylim = c(-maxRank - 0.5, 1),
        clip = "off",
        expand = F
    ) +
    theme(
        legend.position = "bottom",
        panel.border = element_blank()
    ) +
    labs(
        tag = "A"
    )

alleles <- c("MnT", "MT", "FT", "FnT")

hPwDF <- data.frame(
    category = c(trioDF$category, trioDF$category, trioDF$category, trioDF$category),
    snp = c(trioDF$rsid, trioDF$rsid, trioDF$rsid, trioDF$rsid),
    locus = c(trioDF$locus, trioDF$locus, trioDF$locus, trioDF$locus),
    allele = c(rep("MnT", nrow(trioDF)), rep("MT", nrow(trioDF)), rep("FT", nrow(trioDF)), rep("FnT", nrow(trioDF))),
    betaMother = c(trioDF$cmf_bm_moba_pw, trioDF$cmf_bm_moba_pw, trioDF$cmf_bm_moba_pw, trioDF$cmf_bm_moba_pw),
    betaChild = c(trioDF$cmf_bc_moba_pw, trioDF$cmf_bc_moba_pw, trioDF$cmf_bc_moba_pw, trioDF$cmf_bc_moba_pw),
    beta = c(trioDF$h_bmnt_moba_pw, trioDF$h_bmt_moba_pw, trioDF$h_bft_moba_pw, trioDF$h_bfnt_moba_pw),
    se = c(trioDF$h_bmnt_se_moba_pw, trioDF$h_bmt_se_moba_pw, trioDF$h_bft_se_moba_pw, trioDF$h_bfnt_se_moba_pw),
    rank = c(trioDF$rank, trioDF$rank, trioDF$rank, trioDF$rank),
    stringsAsFactors = F
) %>% 
    mutate(
        alleleFactor = factor(allele, levels = alleles),
        scaledBeta = ifelse(category == "MnT_MT", beta / betaMother, beta / betaChild),
        scaledSe = ifelse(category == "MnT_MT", se / betaMother, se / betaChild)
    )


categories <- unique(trioDF$category)

dx <- 0.02

for (i in 1:length(categories)) {
    
    category <- categories[i]
    ranks <- sort(unique(hPwDF$rank[hPwDF$category == category]))
    
    for (j in 1:length(alleles)) {
        
        allele <- alleles[j]
        alelelX <- j /length(alleles) - 0.5/length(alleles)
        
        for (k in 1:length(ranks)) {
            
            rank <- ranks[k]
            
            rankDX <- dx / (length(alleles) + 1) * (k - median(ranks))
            
            x <- alelelX + rankDX
            
            hPwDF$x[hPwDF$category == category & hPwDF$allele == allele & hPwDF$rank == rank] <- x
            
        }
    }
}

plotDF <- hPwDF %>% 
    filter(
        category == "MnT_MT"
    )

maternalPlot <- ggplot() +
    theme_bw() +
    geom_hline(
        yintercept = 0
    ) +
    geom_segment(
        data = plotDF,
        mapping = aes(
            x = x,
            xend = x,
            y = beta - qnorm(0.975) * se + dxCM,
            yend = beta + qnorm(0.975) * se + dxCM
        ),
        alpha = 0.1
    ) +
    geom_point(
        data = plotDF,
        mapping = aes(
            x = x,
            y = beta + dxCM
        ),
        alpha = 0.5
    ) +
    scale_x_continuous(
        limits = c(0, 1),
        expand = c(0, 0),
        breaks = (1:4)/4 - 1/8,
        labels = c("MnT", "MT", "PT", "PnT")
    ) +
    scale_y_continuous(
        name = "Beta [95% CI]"
    ) +
    theme(
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank()
    ) +
    labs(
        tag = "B"
    )

plotDF <- hPwDF %>% 
    filter(
        category == "MT_FT"
    )

fetalPlot <- ggplot() +
    theme_bw() +
    geom_hline(
        yintercept = 0
    ) +
    geom_segment(
        data = plotDF,
        mapping = aes(
            x = x,
            xend = x,
            y = beta - qnorm(0.975) * se + dxCM,
            yend = beta + qnorm(0.975) * se + dxCM
        ),
        alpha = 0.1
    ) +
    geom_point(
        data = plotDF,
        mapping = aes(
            x = x,
            y = beta + dxCM
        ),
        alpha = 0.5
    ) +
    scale_x_continuous(
        breaks = (1:4)/4 - 1/8,
        labels = c("MnT", "MT", "PT", "PnT")
    ) +
    scale_y_continuous(
        name = "Beta [95% CI]"
    ) +
    coord_cartesian(
        xlim = c(0, 1),
        ylim = c(-0.2, 0.3),
        clip = "off",
        expand = F
    ) +
    theme(
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank()
    ) +
    labs(
        tag = "C"
    )


i <- which(trioDF$rsid == "rs2237892")

dx <- 0.2

plotDF <- data.frame(
    x = c(1:4 - dx/2, 1:4 + dx/2),
    beta = c(
        trioDF$h_bmnt_moba_pw[i],
        trioDF$h_bmt_moba_pw[i],
        trioDF$h_bft_moba_pw[i],
        trioDF$h_bfnt_moba_pw[i],
        trioDF$h_bmnt_moba_bw[i],
        trioDF$h_bmt_moba_bw[i],
        trioDF$h_bft_moba_bw[i],
        trioDF$h_bfnt_moba_bw[i]
    ),
    se = c(
        trioDF$h_bmnt_se_moba_pw[i],
        trioDF$h_bmt_se_moba_pw[i],
        trioDF$h_bft_se_moba_pw[i],
        trioDF$h_bfnt_se_moba_pw[i],
        trioDF$h_bmnt_se_moba_bw[i],
        trioDF$h_bmt_se_moba_bw[i],
        trioDF$h_bft_se_moba_bw[i],
        trioDF$h_bfnt_se_moba_bw[i]
    ),
    col = c(
        rep("Placenta Weight (MoBa)", 4),
        rep("Birth Weight (MoBa)", 4)
    ),
    stringsAsFactors = F
)

plotDF$colorFactor <- factor(plotDF$col, levels = c("Placenta Weight (MoBa)", "Birth Weight (MoBa)"))

colors <- c("darkblue", "grey40")

kcnq1Plot <- ggplot() +
    theme_bw() +
    geom_hline(
        yintercept = 0
    ) +
    geom_segment(
        data = plotDF,
        mapping = aes(
            x = x,
            xend = x,
            y = beta - qnorm(0.975) * se,
            yend = beta + qnorm(0.975) * se,
            col = colorFactor
        ),
        alpha = 0.2
    ) +
    geom_segment(
        data = plotDF,
        mapping = aes(
            x = x,
            xend = x,
            y = beta - se,
            yend = beta + se,
            col = colorFactor
        ),
        alpha = 0.2,
        size = 0.8
    ) +
    geom_point(
        data = plotDF,
        mapping = aes(
            x = x,
            y = beta,
            col = colorFactor
        ),
        alpha = 0.5
    ) +
    scale_x_continuous(
        breaks = 1:4,
        labels = c("Maternal\nnon-Transmitted", "Maternal\nTransmitted", "Paternal\nTransmitted", "Paternal\nnon-Transmitted"),
        limits = c(0.5, 4.5),
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        name = "Beta [95% CI]"
    ) +
    scale_color_manual(
        values = colors
    ) + 
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank()
    ) + 
    labs(
        tag = "B"
    )

# Locus Zoom KCNQ1 PW

lzMaxP <- 8
halfWindow <- 150000
ldColors <- brewer.pal(9, "Set1")[c(9, 2, 3, 5, 1)]

locusAnnotationDF <- data.frame(
    name = c("Placenta Weight", "Birth Weight"),
    id = c("rs2237892", "rs234864"),
    pos = c(2839751, 2857297),
    stringsAsFactors = F
)

regionBegin <- min(locusAnnotationDF$pos) - halfWindow
regionEnd <- max(locusAnnotationDF$pos) + halfWindow

recombinationRatesDF <- read.table(
    file = "resources/recombination_rates/genetic_map_GRCh37_chr11.txt.gz",
    header = T,
    sep = "\t",
    stringsAsFactors = F
)

names(recombinationRatesDF) <- c("contig", "position", "rate", "map")

recombinationRatesDF %>%
    filter(
        position >= regionBegin & position <= regionEnd
    ) -> recombinationRatesDF


trioMoba11DF <- read.table(
    file = "/home/marc/Projects/placenta_weight/meta_results/11_pwbw_500k.lm_target.gz",
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>% clean_names() %>% 
    filter(
        contig == "11" & position >= regionBegin & position <= regionEnd
    ) %>% 
    mutate(
        ld = 0
    )

trioMobaPwDF <- trioMoba11DF %>% 
    filter(
        phenotype == "z_placenta_weight"
    )

ld <- fromJSON(file = "resources/ld/targets_IGF2_KCNQ1_ld.json.gz")

variantLD <- ld[["rs2237892"]]

for (i in 1:nrow(trioMobaPwDF)) {
    
    variant2 <- trioMobaPwDF$variant_id[i]
    
    if (variant2 %in% names(variantLD)) {
        
        trioMobaPwDF$ld[i] <- variantLD[[variant2]]
        
    }
}

plotPwDF <- trioMobaPwDF %>% 
    filter(
        !is.na(h_b1_p) & h_b1_p > 0 & !is.na(h_b2_p) & h_b2_p > 0 & !is.na(h_b3_p) & h_b3_p > 0
    ) %>% 
    mutate(
        h_mt = -log10(h_b1_p),
        h_mnt = -log10(h_b2_p),
        h_ft = -log10(h_b3_p)
    ) %>% 
    select(
        variant_id, position, ld, h_mt, h_mnt, h_ft
    ) %>% 
    pivot_longer(
        cols = c("h_mt", "h_mnt", "h_ft"),
        names_to = "haplotype",
        values_to = "logP"
    ) %>% 
    mutate(
        haplotypeFactor = factor(haplotype, levels = c("h_mnt", "h_mt", "h_ft")),
        ldFactor = case_when(
            ld < 0.2 ~ "[0.0  0.2]",
            ld >= 0.2 & ld < 0.4 ~ "[0.2  0.4]",
            ld >= 0.4 & ld < 0.6 ~ "[0.4  0.6]",
            ld >= 0.6 & ld < 0.8 ~ "[0.6  0.8]",
            ld >= 0.8 ~ "[0.8  1.0]"
        ),
        ldFactor = factor(ldFactor, c("[0.0  0.2]", "[0.2  0.4]", "[0.4  0.6]", "[0.6  0.8]", "[0.8  1.0]"))
    ) %>% 
    arrange(
        ldFactor
    )

levels(plotPwDF$haplotypeFactor) <- c("MnT", "MT", "PT")

annotationPwDF1 <- plotPwDF %>% 
    filter(
        variant_id == locusAnnotationDF$id[1]
    )

annotationPwDF2 <- plotPwDF %>% 
    filter(
        variant_id == locusAnnotationDF$id[2]
    )

recombinationRatesDF <- recombinationRatesDF %>% 
    mutate(
        rateScaled = 0.8 * rate / max(rate) * lzMaxP
    )

labelPwDF <- locusAnnotationDF %>% 
    left_join(
        plotPwDF %>% 
            filter(
                haplotype == "h_mt"
            ),
        by = c(id = "variant_id")
    ) %>% 
    mutate(
        dx = case_when(
            id == "rs2237892" ~ - 50000,
            id == "rs234864" ~ + 50000
        ),
        y = 7,
        hjust = case_when(
            id == "rs2237892" ~ 1,
            id == "rs234864" ~ 0
        )
    )


kcnq1LzPlot <- ggplot() +
    theme_bw() + 
    geom_line(
        data = recombinationRatesDF,
        mapping = aes(
            x = position,
            y = rateScaled
        ),
        col = "blue3"
    ) +
    geom_segment(
        data = labelPwDF,
        mapping = aes(
            x = position,
            xend = position + dx,
            y = logP,
            yend = y
        ),
        linetype = "dotted"
    ) +
    geom_label(
        data = labelPwDF,
        mapping = aes(
            x = position + dx,
            y = y,
            nudge_x = dx,
            label = name,
            hjust = hjust
        ),
        vjust = 0.5
    ) +
    geom_point(
        data = plotPwDF,
        mapping = aes(
            x = position,
            y = logP,
            col = ldFactor
        ),
        size = 2,
        alpha = 0.8
    ) +
    geom_point(
        data = annotationPwDF2,
        mapping = aes(
            x = position,
            y = logP
        ),
        size = 3,
        color = "black"
    ) +
    geom_point(
        data = annotationPwDF2,
        mapping = aes(
            x = position,
            y = logP
        ),
        size = 2,
        color = ldColors[1]
    ) +
    geom_point(
        data = annotationPwDF1,
        mapping = aes(
            x = position,
            y = logP
        ),
        shape = 18,
        col = "grey20",
        size = 5
    ) +
    geom_point(
        data = annotationPwDF1,
        mapping = aes(
            x = position,
            y = logP
        ),
        shape = 18,
        col = "red3",
        size = 4
    ) +
    scale_color_manual(
        name = "R2",
        values = ldColors,
        guide = guide_legend(
            override.aes = list(
                size = 4
            )
        )
    ) + 
    scale_x_continuous(
        breaks = regionBegin + halfWindow,
        labels = glue("{round((regionBegin + halfWindow)/1000000, 2)} Mb")
    ) + 
    scale_y_continuous(
        name = "Placenta Weight [MoBa]",
        limits = c(0, lzMaxP)
    ) +
    coord_cartesian(
        clip = "off"
    ) +
    facet_grid(
        . ~ haplotypeFactor
    ) +
    theme(
        legend.position = "top",
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(
            fill = "white",
            colour = NA
        )
    ) + 
    labs(
        tag = "C"
    )


# Locus Zoom KCNQ1 BW

locusAnnotationDF <- data.frame(
    name = c("Placenta Weight", "Birth Weight"),
    id = c("rs2237892", "rs234864"),
    pos = c(2839751, 2857297),
    stringsAsFactors = F
)

regionBegin <- min(locusAnnotationDF$pos) - halfWindow
regionEnd <- max(locusAnnotationDF$pos) + halfWindow

recombinationRatesDF <- read.table(
    file = "resources/recombination_rates/genetic_map_GRCh37_chr11.txt.gz",
    header = T,
    sep = "\t",
    stringsAsFactors = F
)

names(recombinationRatesDF) <- c("contig", "position", "rate", "map")

recombinationRatesDF %>%
    filter(
        position >= regionBegin & position <= regionEnd
    ) -> recombinationRatesDF


trioMobaBwDF <- trioMoba11DF %>% 
    filter(
        phenotype == "z_weight0"
    )

ld <- fromJSON(file = "resources/ld/targets_IGF2_KCNQ1_ld.json.gz")

variantLD <- ld[["rs2237892"]]

for (i in 1:nrow(trioMobaBwDF)) {
    
    variant2 <- trioMobaBwDF$variant_id[i]
    
    if (variant2 %in% names(variantLD)) {
        
        trioMobaBwDF$ld[i] <- variantLD[[variant2]]
        
    }
}

plotBwDF <- trioMobaBwDF %>% 
    filter(
        !is.na(h_b1_p) & h_b1_p > 0 & !is.na(h_b2_p) & h_b2_p > 0 & !is.na(h_b3_p) & h_b3_p > 0
    ) %>% 
    mutate(
        h_mt = -log10(h_b1_p),
        h_mnt = -log10(h_b2_p),
        h_ft = -log10(h_b3_p)
    ) %>% 
    select(
        variant_id, position, ld, h_mt, h_mnt, h_ft
    ) %>% 
    pivot_longer(
        cols = c("h_mt", "h_mnt", "h_ft"),
        names_to = "haplotype",
        values_to = "logP"
    ) %>% 
    mutate(
        haplotypeFactor = factor(haplotype, levels = c("h_mnt", "h_mt", "h_ft")),
        ldFactor = case_when(
            ld < 0.2 ~ "[0.0  0.2]",
            ld >= 0.2 & ld < 0.4 ~ "[0.2  0.4]",
            ld >= 0.4 & ld < 0.6 ~ "[0.4  0.6]",
            ld >= 0.6 & ld < 0.8 ~ "[0.6  0.8]",
            ld >= 0.8 ~ "[0.8  1.0]"
        ),
        ldFactor = factor(ldFactor, c("[0.0  0.2]", "[0.2  0.4]", "[0.4  0.6]", "[0.6  0.8]", "[0.8  1.0]"))
    ) %>% 
    arrange(
        ldFactor
    )

levels(plotBwDF$haplotypeFactor) <- c("MnT", "MT", "PT")

annotationBwDF1 <- plotBwDF %>% 
    filter(
        variant_id == locusAnnotationDF$id[1]
    )

annotationBwDF2 <- plotBwDF %>% 
    filter(
        variant_id == locusAnnotationDF$id[2]
    )

recombinationRatesDF <- recombinationRatesDF %>% 
    mutate(
        rateScaled = 0.8 * rate / max(rate) * lzMaxP
    )

labelBwDF <- locusAnnotationDF %>% 
    left_join(
        plotBwDF %>% 
            filter(
                haplotype == "h_mt"
            ),
        by = c(id = "variant_id")
    ) %>% 
    mutate(
        dx = case_when(
            id == "rs2237892" ~ - 50000,
            id == "rs234864" ~ + 50000
        ),
        y = 5,
        hjust = case_when(
            id == "rs2237892" ~ 1,
            id == "rs234864" ~ 0
        )
    )


bwLzPlot <- ggplot() +
    theme_bw() + 
    geom_line(
        data = recombinationRatesDF,
        mapping = aes(
            x = position,
            y = rateScaled
        ),
        col = "blue3"
    ) +
    geom_segment(
        data = labelBwDF,
        mapping = aes(
            x = position,
            xend = position + dx,
            y = logP,
            yend = y
        ),
        linetype = "dotted"
    ) +
    geom_label(
        data = labelBwDF,
        mapping = aes(
            x = position + dx,
            y = y,
            nudge_x = dx,
            label = name,
            hjust = hjust
        ),
        vjust = 0.5
    ) +
    geom_point(
        data = plotBwDF,
        mapping = aes(
            x = position,
            y = logP,
            col = ldFactor
        ),
        size = 2,
        alpha = 0.8
    ) +
    geom_point(
        data = annotationBwDF2,
        mapping = aes(
            x = position,
            y = logP
        ),
        size = 3,
        color = "black"
    ) +
    geom_point(
        data = annotationBwDF2,
        mapping = aes(
            x = position,
            y = logP
        ),
        size = 2,
        color = ldColors[1]
    ) +
    geom_point(
        data = annotationBwDF1,
        mapping = aes(
            x = position,
            y = logP
        ),
        shape = 18,
        col = "grey20",
        size = 5
    ) +
    geom_point(
        data = annotationBwDF1,
        mapping = aes(
            x = position,
            y = logP
        ),
        shape = 18,
        col = "red3",
        size = 4
    ) +
    scale_color_manual(
        name = "R2",
        values = ldColors,
        guide = guide_legend(
            override.aes = list(
                size = 4
            )
        )
    ) + 
    scale_x_continuous(
        breaks = regionBegin + halfWindow,
        labels = glue("{round((regionBegin + halfWindow)/1000000, 2)} Mb")
    ) + 
    scale_y_continuous(
        name = "Birth Weight [MoBa]",
        limits = c(0, lzMaxP)
    ) +
    coord_cartesian(
        clip = "off"
    ) +
    facet_grid(
        . ~ haplotypeFactor
    ) +
    theme(
        legend.position = "top",
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(
            fill = "white",
            colour = NA
        )
    ) + 
    labs(
        tag = "D"
    )


pdf(
    file = glue("docs/figures/trio/fig_4_allele_contribution_15_09_22.pdf"),
    width = unit(12, "cm"),
    height = unit(24.7 / 3, "cm"), 
    pointsize = 12
)
(heatMap | (kcnq1Plot / kcnq1LzPlot / bwLzPlot)) + 
    plot_layout(
        widths = c(2, 1.2)
    )
dummy <- dev.off()

