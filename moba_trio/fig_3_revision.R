
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
library(here)


# Namespace conflits

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")


# Load data

trioDF <- read.table(
    file = file.path(here(), "moba_trio/resources/merge_decode_06.07.22.txt"),
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>% 
    clean_names()

metaResultsDF <- read.table(
    file = file.path(here(), "moba_trio/resources/meta_results_07.02.23.txt"),
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>% 
    clean_names() %>% 
    mutate(
        ea = toupper(ea),
        oa = toupper(oa)
    )


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
    
    ea <- metaResultsDF$ea[metaResultsDF$rsid == rsid]
    
    if (length(ea) != 1) {
        
        stop(glue("Effect allele for {rsid} not found."))
        
    }
    
    oa <- metaResultsDF$oa[metaResultsDF$rsid == rsid]
    
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


# Merge with meta

rows_before <- nrow(trioDF)

trioDF <- trioDF %>% 
    left_join(
        metaResultsDF,
        by = "rsid"
    )

rows_after <- nrow(trioDF)

if (sum(trioDF$oa != trioDF$other_allele) > 0) {
    
    stop("Other allele mismatch")
    
}
if (sum(trioDF$ea != trioDF$effect_allele) > 0) {
    
    stop("Effect allele mismatch")
    
}
if (rows_before != rows_after) {
    
    stop("Duplicates introduced when merging")
    
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


# Effect size and se of all variants

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
                rep("Placental Weight (MoBa)", 4),
                rep("Birth Weight (MoBa)", 4),
                rep("Birth Weight (deCode)", 3)
            ),
            stringsAsFactors = F
        ) %>% 
            mutate(
                ymin = beta - qnorm(0.975) * se,
                ymax = beta + qnorm(0.975) * se
            )
        
        plotDF$colorFactor <- factor(plotDF$col, levels = c("Placental Weight (MoBa)", "Birth Weight (MoBa)", "Birth Weight (deCode)"))
        
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
                rep("Placental Weight (MoBa)", 4),
                rep("Birth Weight (MoBa)", 4)
            ),
            stringsAsFactors = F
        ) %>% 
            mutate(
                ymin = beta - qnorm(0.975) * se,
                ymax = beta + qnorm(0.975) * se
            )
        
        plotDF$colorFactor <- factor(plotDF$col, levels = c("Placental Weight (MoBa)", "Birth Weight (MoBa)"))
        
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
            rep("Placental Weight (MoBa)", 4),
            rep("Birth Weight (MoBa)", 4)
        ),
        stringsAsFactors = F
    ) %>% 
        mutate(
            ymin = beta - qnorm(0.975) * se,
            ymax = beta + qnorm(0.975) * se
        )
    
    plotDF$colorFactor <- factor(plotDF$col, levels = c("Placental Weight (MoBa)", "Birth Weight (MoBa)"))
    
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
    file = glue("{here()}/moba_trio/figures/haplotypes_supplementary_03.02.22.pdf"),
    width = unit(21, "cm"),
    height = unit(3 * 29.7, "cm"), 
    pointsize = 10
)
grid.draw(plotPatchwork)
dummy <- dev.off()

pdf(
    file = glue("{here()}/moba_trio/figures/haplotypes_supplementary_1_03.02.22.pdf"),
    width = unit(21, "cm"),
    height = unit(29.7, "cm"), 
    pointsize = 10
)
grid.draw(plotPatchwork1)
dummy <- dev.off()

pdf(
    file = glue("{here()}/moba_trio/figures/haplotypes_supplementary_2_03.02.22.pdf"),
    width = unit(21, "cm"),
    height = unit(29.7, "cm"), 
    pointsize = 10
)
grid.draw(plotPatchwork2)
dummy <- dev.off()
pdf(
    file = glue("{here()}/moba_trio/figures/haplotypes_supplementary_3_03.02.22.pdf"),
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
        maternal_fetal_classification_factor, pw_bw_classification, rank
    ) %>% 
    mutate(
        rank = row_number() + 0.5 * (as.numeric(maternal_fetal_classification_factor) - 1)
    ) %>% 
    separate(
        locus, 
        into = c("locus"), 
        extra = "drop"
    )

m_min <- min(trioDF$rank[trioDF$maternal_fetal_classification == "Maternal"])
m_max <- max(trioDF$rank[trioDF$maternal_fetal_classification == "Maternal"])
fm_min <- min(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal & Maternal"])
fm_max <- max(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal & Maternal"])
fmo_min <- min(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal & Maternal opposite directions"])
fmo_max <- max(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal & Maternal opposite directions"])
f_pw_min <- min(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal" & trioDF$pw_bw_classification == "Predominantly or only PW"])
f_pw_max <- max(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal" & trioDF$pw_bw_classification == "Predominantly or only PW"])
f_pwbw_min <- min(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal" & trioDF$pw_bw_classification == "PW & BW same direction"])
f_pwbw_max <- max(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal" & trioDF$pw_bw_classification == "PW & BW same direction"])
f_pwbwo_min <- min(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal" & trioDF$pw_bw_classification == "PW & BW opposite directions"])
f_pwbwo_max <- max(trioDF$rank[trioDF$maternal_fetal_classification == "Fetal" & trioDF$pw_bw_classification == "PW & BW opposite directions"])
u_pw_min <- min(trioDF$rank[trioDF$maternal_fetal_classification == "Unclassified" & trioDF$pw_bw_classification == "Predominantly or only PW"])
u_pw_max <- max(trioDF$rank[trioDF$maternal_fetal_classification == "Unclassified" & trioDF$pw_bw_classification == "Predominantly or only PW"])
u_pwbw_min <- min(trioDF$rank[trioDF$maternal_fetal_classification == "Unclassified" & trioDF$pw_bw_classification == "PW & BW same direction"])
u_pwbw_max <- max(trioDF$rank[trioDF$maternal_fetal_classification == "Unclassified" & trioDF$pw_bw_classification == "PW & BW same direction"])


separator_top <- c(m_min, fm_min, fmo_min, f_pw_min, f_pwbw_min, f_pwbwo_min, u_pw_min, u_pwbw_min)
separator_top <- -separator_top + 0.5

separator_bottom <- c(m_max, fm_max, fmo_max, f_pw_max, f_pwbw_max, f_pwbwo_max, u_pw_max, u_pwbw_max)
separator_bottom <- -separator_bottom - 0.5

separators <- c(separator_top, separator_bottom)

maxRank <- max(trioDF$rank)

fill_variable <- "beta" # Use "t" or "beta"

betaDF <- data.frame(
    rank = c(trioDF$rank, trioDF$rank, trioDF$rank, trioDF$rank, 
             trioDF$rank, trioDF$rank, trioDF$rank,
             trioDF$rank, trioDF$rank, trioDF$rank, trioDF$rank),
    variable = c(rep("mnt_moba_pw", nrow(trioDF)), rep("mt_moba_pw", nrow(trioDF)), rep("ft_moba_pw", nrow(trioDF)), rep("fnt_moba_pw", nrow(trioDF)), 
                 rep("mnt_decode_bw", nrow(trioDF)), rep("mt_decode_bw", nrow(trioDF)), rep("ft_decode_bw", nrow(trioDF)),
                 rep("f_meta_pw", nrow(trioDF)), rep("m_meta_pw", nrow(trioDF)), rep("f_meta_bw", nrow(trioDF)), rep("m_meta_bw", nrow(trioDF))),
    beta = c(trioDF$h_bmnt_moba_pw, trioDF$h_bmt_moba_pw, trioDF$h_bft_moba_pw, trioDF$h_bfnt_moba_pw, 
             trioDF$decode_mnt_beta, trioDF$decode_mt_beta, trioDF$decode_ft_beta,
             trioDF$meta_beta_fetal_pw, trioDF$meta_beta_maternal_pw, trioDF$meta_beta_fetal_bw, trioDF$meta_beta_maternal_bw),
    se = c(trioDF$h_bmnt_se_moba_pw, trioDF$h_bmt_se_moba_pw, trioDF$h_bft_se_moba_pw, trioDF$h_bfnt_se_moba_pw, 
           trioDF$decode_mnt_se, trioDF$decode_mt_se, trioDF$decode_ft_se,
           trioDF$meta_se_fetal_pw, trioDF$meta_se_maternal_pw, trioDF$meta_se_fetal_bw, trioDF$meta_se_maternal_bw),
    p = c(trioDF$h_bmnt_p_moba_pw, trioDF$h_bmt_p_moba_pw, trioDF$h_bft_p_moba_pw, trioDF$h_bfnt_p_moba_pw, 
          trioDF$decode_mnt_p, trioDF$decode_mt_p, trioDF$decode_ft_p,
          trioDF$meta_p_fetal_pw, trioDF$meta_p_maternal_pw, trioDF$meta_p_fetal_bw, trioDF$meta_p_maternal_bw),
    stringsAsFactors = F
) %>% 
    mutate(
        variableFactor = factor(variable, levels = c("m_meta_pw", "f_meta_pw", "mnt_moba_pw", "mt_moba_pw", "ft_moba_pw", "fnt_moba_pw", "mnt_decode_bw", "mt_decode_bw", "ft_decode_bw", "m_meta_bw", "f_meta_bw")),
        x = as.numeric(variableFactor),
        y = -rank,
        t = beta / se,
        logP = -log10(p),
        scale = case_when(
            logP > -log10(5e-8) ~ 1,
            logP > 6 ~ 0.8,
            logP > 3 ~ 0.6,
            logP > -log10(0.05) ~ 0.4,
            T ~ 0.2
        )
    )

missingDF <- betaDF %>% 
    filter(
        is.na(beta)
    )

betaDF <- betaDF %>% 
    filter(
        !is.na(beta)
    )

if (fill_variable == "beta") {
    
    fill_limits <- max(abs(betaDF$beta[!is.na(betaDF$beta)]))
    betaDF$fill <- betaDF$beta
    legend_title <- "Beta"
    
} else {
    
    fill_limits <- max(abs(betaDF$t[!is.na(betaDF$t)]))
    betaDF$fill <- betaDF$t
    legend_title <- "Beta / se"
    
}

first_box_x <- 6.5
second_box_x <- 9.5

legendDF <- data.frame(
    x = c(rep(first_box_x, 3), rep(second_box_x, 3)),
    y = -u_pwbw_max - 1.5 - c(0:2, 0:2),
    label = c("<i>P</i> < 5 × 10<sup>-8</sup>", "<i>P</i> < 10<sup>-6</sup>", "<i>P</i> < 10<sup>-3</sup>", "<i>P</i> < 0.05", "<i>P</i> >= 0.05", "No proxy"),
    scale = c(1, 0.8, 0.6, 0.4, 0.2, 0)
)


column_labels <- c("M", "F", "MnT", "MT", "PT", "PnT", "MnT", "MT", "PT", "M", "F")
n_columns <- length(column_labels)

heatMap <- ggplot() +
    theme_void() +
    geom_text(
        mapping = aes(
            x = 1:n_columns,
            y = 0,
            label = column_labels
        ),
        size = 2.5,
        col = "grey40",
        hjust = 0.5,
        vjust = 0
    ) +
    geom_richtext(
        data = trioDF,
        mapping = aes(
            x = 0.4,
            y = -rank,
            label = glue("<i>{locus}</i> ({rsid})")
        ),
        size = 2.5,
        col = "grey30",
        hjust = 1,
        vjust = 0.5,
        fill = NA, 
        label.color = NA
    ) +
    geom_segment(
        data = missingDF,
        mapping = aes(
            x = x - 0.15,
            xend = x + 0.15,
            y = y - 0.15,
            yend = y + 0.15
        ),
        col = "grey70"
    ) +
    geom_segment(
        data = missingDF,
        mapping = aes(
            x = x - 0.15,
            xend = x + 0.15,
            y = y + 0.15,
            yend = y - 0.15
        ),
        col = "grey70"
    ) +
    geom_tile(
        data = betaDF,
        mapping = aes(
            x = x,
            y = y,
            fill = fill,
            width = scale,
            height = scale
        )
    ) +
    geom_segment(
        mapping = aes(
            x = c(0.5, 2.5, 6.5, 9.5, n_columns + 0.5),
            xend = c(0.5, 2.5, 6.5, 9.5, n_columns + 0.5),
            y = -0.25,
            yend = -m_max - 0.5
        ),
        col = "grey70"
    ) +
    geom_segment(
        mapping = aes(
            x = c(0.5, 2.5, 6.5, 9.5, n_columns + 0.5),
            xend = c(0.5, 2.5, 6.5, 9.5, n_columns + 0.5),
            y = -fm_min + 0.5,
            yend = -fm_max - 0.5
        ),
        col = "grey70"
    ) +
    geom_segment(
        mapping = aes(
            x = c(0.5, 2.5, 6.5, 9.5, n_columns + 0.5),
            xend = c(0.5, 2.5, 6.5, 9.5, n_columns + 0.5),
            y = -fmo_min + 0.5,
            yend = -fmo_max - 0.5
        ),
        col = "grey70"
    ) +
    geom_segment(
        mapping = aes(
            x = c(0.5, 2.5, 6.5, 9.5, n_columns + 0.5),
            xend = c(0.5, 2.5, 6.5, 9.5, n_columns + 0.5),
            y = -f_pw_min + 0.5,
            yend = -f_pwbwo_max - 0.5
        ),
        col = "grey70"
    ) +
    geom_segment(
        mapping = aes(
            x = c(0.5, 2.5, 6.5, 9.5, n_columns + 0.5),
            xend = c(0.5, 2.5, 6.5, 9.5, n_columns + 0.5),
            y = -u_pw_min + 0.5,
            yend = -u_pwbw_max - 0.5
        ),
        col = "grey60"
    ) +
    geom_segment(
        mapping = aes(
            x = 0.5,
            xend = n_columns + 0.6,
            y = separators,
            yend = separators
        ),
        col = "grey60"
    ) +
    geom_text(
        mapping = aes(
            x = -3.4,
            y = c(-(m_min + m_max) / 2, -(fm_min + fmo_max) / 2, -(f_pw_min + f_pwbwo_max) / 2, -(u_pw_min + u_pwbw_max) / 2),
            label = c("Maternal", "F & M", "Fetal", "Unclassified")
        ),
        size = 3.5,
        col = "grey30",
        hjust = 0.5,
        vjust = 1,
        angle = 90
    ) +
    geom_text(
        mapping = aes(
            x = n_columns + 0.75,
            y = c(-(m_min + m_max) / 2, -(fm_min + fm_max) / 2, -(fmo_min + fmo_max) / 2, -(f_pw_min + f_pw_max) / 2, -(f_pwbw_min + f_pwbw_max) / 2, -(f_pwbwo_min + f_pwbwo_max) / 2, -(u_pw_min + u_pw_max) / 2, -(u_pwbw_min + u_pwbw_max) / 2),
            label = c("PW & BW", "Predominantly PW", "PW & BW", "Predominantly PW", "PW & BW", "PW & BW opposite", "Predominantly PW", "PW & BW")
        ),
        size = 2.5,
        col = "grey30",
        hjust = 0,
        vjust = 0.5
    ) +
    geom_richtext(
        mapping = aes(
            x = c(1.5, 4.5, 8, 10.5),
            y = 1.2,
            label = c("PW<br>Meta", "PW<br>Transmission", "BW<br>Transmission<sup>1</sup>", "BW<br>Meta<sup>2</sup>")
        ),
        size = 3.2,
        col = "grey30",
        hjust = 0.5,
        vjust = 0,
        fill = NA, 
        label.color = NA
    ) +
    geom_richtext(
        mapping = aes(
            x = -3.4,
            y = c(-u_pwbw_max - 2.5, -u_pwbw_max - 3.5),
            label = c("<sup>1</sup> Juliusdottir et al.", "<sup>2</sup> Warrington et al.")
        ),
        size = 3.2,
        col = "grey30",
        hjust = 0,
        vjust = 0,
        fill = NA, 
        label.color = NA
    ) +
    geom_segment(
        data = missingDF,
        mapping = aes(
            x = second_box_x - 0.15,
            xend = second_box_x + 0.15,
            y = -u_pwbw_max - 3.5 - 0.15,
            yend = -u_pwbw_max - 3.5 + 0.15
        ),
        col = "grey70"
    ) +
    geom_segment(
        data = missingDF,
        mapping = aes(
            x = second_box_x - 0.15,
            xend = second_box_x + 0.15,
            y = -u_pwbw_max - 3.5 + 0.15,
            yend = -u_pwbw_max - 3.5 - 0.15
        ),
        col = "grey70"
    ) +
    geom_tile(
        data = legendDF,
        mapping = aes(
            x = x,
            y = y,
            width = scale,
            height = scale
        ),
        fill = "grey80"
    ) +
    geom_richtext(
        data = legendDF,
        mapping = aes(
            x = x + 0.52,
            y = y,
            label = label
        ),
        hjust = 0,
        size = 3.2,
        fill = NA, 
        label.color = NA
    ) +
    scale_fill_scico(
        name = legend_title,
        palette = "vik",
        limits = c(-fill_limits, fill_limits)
    ) +
    coord_cartesian(
        xlim = c(-4, 14.8),
        ylim = c(-maxRank - 0.5 - 4, 0),
        clip = "off",
        expand = F
    ) +
    theme(
        legend.position = c(.35, 0.033), 
        legend.direction = "horizontal",
        panel.border = element_blank()
    ) +
    labs(
        tag = "a"
    )

# heatMap

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
        tag = "b"
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
        tag = "c"
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
        rep("Placental Weight (MoBa)", 4),
        rep("Birth Weight (MoBa)", 4)
    ),
    stringsAsFactors = F
)

plotDF$colorFactor <- factor(plotDF$col, levels = c("Placental Weight (MoBa)", "Birth Weight (MoBa)"))

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
        linewidth = 0.8
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
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()
    ) + 
    labs(
        tag = "b"
    )


# Gene coordinates

regionBegin <- min(locusAnnotationDF$pos) - halfWindow
regionEnd <- max(locusAnnotationDF$pos) + halfWindow

minWidth <- (regionEnd - regionBegin)/40

geneCoordinatesDF <- data.frame(
  name = c("KCNQ1", "KCNQ1OT1", "COX6CP18", "KCNQ1-AS1", "KCNQ1DN", "CDKN1C"),
  start = c(2465914, 2629558, 2792586, 2861365,	2891263, 2904443),
  end = c(2870339, 2721224,	2792807, 2882798, 2893335, 2904443)
) %>% 
  filter(
    name %in% c("KCNQ1", "KCNQ1-AS1", "CDKN1C")
  )

geneCoordinatesDF <- geneCoordinatesDF %>% 
  filter(
    end > regionBegin & start < regionEnd
  ) %>% 
  mutate(
    start = ifelse(start < regionBegin, regionBegin, start),
    end = ifelse(end < regionBegin, regionBegin, end)
  )

geneCoordinatesPlotDF <- geneCoordinatesDF %>% 
  mutate(
    names_factor = factor(name, levels = geneCoordinatesDF$name),
    y = (nrow(geneCoordinatesDF) - as.numeric(names_factor)) / nrow(geneCoordinatesDF) - 1,
    end = ifelse(end - start < minWidth, start + minWidth, end)
  )

geneCoordinatesPlotDF <- rbind(
  geneCoordinatesPlotDF %>% 
    mutate(
      allele_transmission = "MnT"
    ),
  geneCoordinatesPlotDF %>% 
    mutate(
      allele_transmission = "MT"
    ),
  geneCoordinatesPlotDF %>% 
    mutate(
      allele_transmission = "PT"
    )
)

exonDF <- rbind(
  data.frame(
    name = "KCNQ1",
    start = c(2466221, 2549158,	2591858,	2592555,	2593243,	2594076,	2604665,	2606442,	2608800,	2609943,	2683191,	2790074,	2797190,	2798216,	2799206,	2868997),
    end = c(2466714,	2549248,	2591984,	2592633,	2593339,	2594216,	2604775,	2606537,	2608922,	2610084,	2683311,	2790149,	2797284,	2798262,	2799267,	2870339)
  ),
  data.frame(
    name = "KCNQ1-AS1",
    start = c(2882798,	2880169,	2861672),
    end = c(2882220	,2879973,	2861365)
  ),
  data.frame(
    name = "CDKN1C",
    start = c(2907111,	2905364,	2905145),
    end = c(2905900,	2905229,	2904443
    )
  )
)

exonDF <- exonDF %>% 
  filter(
    end > regionBegin & start < regionEnd
  ) %>% 
  mutate(
    start = ifelse(start < regionBegin, regionBegin, start),
    end = ifelse(end < regionBegin, regionBegin, end)
  )

exonPlotDF <- exonDF %>% 
  mutate(
    names_factor = factor(name, levels = geneCoordinatesDF$name),
    y = (nrow(geneCoordinatesDF) - as.numeric(names_factor)) / nrow(geneCoordinatesDF) - 1,
    end = ifelse(end - start < minWidth, start + minWidth, end)
  )

exonPlotDF <- rbind(
  exonPlotDF %>% 
    mutate(
      allele_transmission = "MnT"
    ),
  exonPlotDF %>% 
    mutate(
      allele_transmission = "MT"
    ),
  exonPlotDF %>% 
    mutate(
      allele_transmission = "PT"
    )
)

# Locus Zoom KCNQ1 PW

lzMaxP <- 7
halfWindow <- 150000
ldColors <- brewer.pal(9, "Set1")[c(9, 2, 3, 5, 1)]

locusAnnotationDF <- data.frame(
    name = c("Placental Weight", "Birth Weight"),
    id = c("rs2237892", "rs234864"),
    pos = c(2839751, 2857297),
    stringsAsFactors = F
)

regionBegin <- min(locusAnnotationDF$pos) - halfWindow
regionEnd <- max(locusAnnotationDF$pos) + halfWindow

recombinationRatesDF <- read.table(
    file = glue("{here()}/moba_trio/resources/genetic_map_GRCh37_chr11.txt.gz"),
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
    file = glue("{here()}/moba_trio/resources/kcnq1_lz.gz"),
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>% 
    mutate(
        ld = 0
    )

trioMobaPwDF <- trioMoba11DF %>% 
    filter(
        phenotype == "z_placenta_weight"
    )

ld <- fromJSON(file = glue("{here()}/moba_trio/resources/targets_IGF2_KCNQ1_ld.json.gz"))

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

levels(plotPwDF$haplotypeFactor) <- c("Maternal\nnon-Transmitted", "Maternal\nTransmitted", "Paternal\nTransmitted")

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
        vjust = 0.5,
        size = 2.5
    ) +
    geom_point(
        data = plotPwDF,
        mapping = aes(
            x = position,
            y = logP,
            col = ldFactor
        ),
        size = 1,
        alpha = 0.8
    ) +
    geom_point(
        data = annotationPwDF2,
        mapping = aes(
            x = position,
            y = logP
        ),
        size = 2,
        color = "black"
    ) +
    geom_point(
        data = annotationPwDF2,
        mapping = aes(
            x = position,
            y = logP
        ),
        size = 1,
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
        size = 4
    ) +
    geom_point(
        data = annotationPwDF1,
        mapping = aes(
            x = position,
            y = logP
        ),
        shape = 18,
        col = "red3",
        size = 3
    ) +
  geom_segment(
    data = geneCoordinatesPlotDF,
    mapping = aes(
      x = start,
      xend = end,
      y = y,
      yend = y
    ),
    col = "grey30",
    size = 0.3
  )  +
  geom_rect(
    data = exonPlotDF,
    mapping = aes(
      xmin = start,
      xmax = end,
      ymin = y - 1/20,
      ymax = y + 1/20
    ),
    col = "grey30",
    fill = "grey90",
    size = 0.3
  ) +
  geom_richtext(
    data = geneCoordinatesPlotDF,
    mapping = aes(
      x = end + minWidth/2,
      y = y,
      label = glue("<i>{name}</i>")
    ),
    size = 2,
    col = "grey30",
    hjust = 0,
    vjust = 0.5,
    fill = NA, 
    label.color = NA
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
        name = "Placental Weight [MoBa]",
        limits = c(-1.05, lzMaxP),
        expand = c(0, 0)
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
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text = element_blank()
    ) + 
    labs(
        tag = "c"
    )


# Locus Zoom KCNQ1 BW

locusAnnotationDF <- data.frame(
    name = c("Placental Weight", "Birth Weight"),
    id = c("rs2237892", "rs234864"),
    pos = c(2839751, 2857297),
    stringsAsFactors = F
)

recombinationRatesDF <- read.table(
  file = glue("{here()}/moba_trio/resources/genetic_map_GRCh37_chr11.txt.gz"),
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

ld <- fromJSON(file = glue("{here()}/moba_trio/resources/targets_IGF2_KCNQ1_ld.json.gz"))

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

levels(plotBwDF$haplotypeFactor) <- c("Maternal\nnon-Transmitted", "Maternal\nTransmitted", "Paternal\nTransmitted")

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
        vjust = 0.5,
        size = 2.5
    ) +
    geom_point(
        data = plotBwDF,
        mapping = aes(
            x = position,
            y = logP,
            col = ldFactor
        ),
        size = 1,
        alpha = 0.8
    ) +
    geom_point(
        data = annotationBwDF2,
        mapping = aes(
            x = position,
            y = logP
        ),
        size = 2,
        color = "black"
    ) +
    geom_point(
        data = annotationBwDF2,
        mapping = aes(
            x = position,
            y = logP
        ),
        size = 1,
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
        size = 4
    ) +
    geom_point(
        data = annotationBwDF1,
        mapping = aes(
            x = position,
            y = logP
        ),
        shape = 18,
        col = "red3",
        size = 3
    ) +
  geom_segment(
    data = geneCoordinatesPlotDF,
    mapping = aes(
      x = start,
      xend = end,
      y = y,
      yend = y
    ),
    col = "grey30",
    size = 0.3
  )  +
  geom_rect(
    data = exonPlotDF,
    mapping = aes(
      xmin = start,
      xmax = end,
      ymin = y - 1/20,
      ymax = y + 1/20
    ),
    col = "grey30",
    fill = "grey90",
    size = 0.3
  ) +
  geom_richtext(
    data = geneCoordinatesPlotDF,
    mapping = aes(
      x = end + minWidth/2,
      y = y,
      label = glue("<i>{name}</i>")
    ),
    size = 2,
    col = "grey30",
    hjust = 0,
    vjust = 0.5,
    fill = NA, 
    label.color = NA
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
        limits = c(-1.05, lzMaxP),
        expand = c(0, 0)
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
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(
          fill = "white",
          colour = NA
        )
    ) + 
    labs(
        tag = "d"
    )
  



pdf(
    file = glue("{here()}/moba_trio/figures//fig_3_mode_of_association_06_07_23_{fill_variable}_scaled.pdf"),
    width = unit(12, "cm"),
    height = unit(24.7 / 3, "cm"), 
    pointsize = 12
)
(heatMap | ((kcnq1Plot / ((kcnq1LzPlot / bwLzPlot) + plot_layout(guides = "collect") & theme(legend.position = 'top')))) + plot_layout(heights = c(0.7, 2.3))) + 
    plot_layout(
        widths = c(2, 1.2),
        guides = "keep"
    )
dummy <- dev.off()

