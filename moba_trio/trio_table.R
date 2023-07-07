
##
#
# This script gets the data for the tables of the trio and haplotype analyses.
#
##


# Libraries

library(conflicted)
library(glue)
library(janitor)
library(tidyr)
library(dplyr)
library(rjson)
library(curl)
library(jsonlite)
library(httr)
library(stringr)
library(glue)
library(LDlinkR)


# Namespace conflits

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("fromJSON", "rjson")
conflict_prefer("toJSON", "rjson")


# Ensembl API

server <- "https://grch37.rest.ensembl.org"


getGeneIdentifier <- function(
    geneName
) {
    
    ext <- paste0("/xrefs/symbol/homo_sapiens/", geneName, "?")
    
    geneJson <- NULL
    
    for (try in 1:100) {
        
        tryCatch(
            {
                r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
                
                responseCode <- status_code(r)
                
                if (responseCode == 200) {
                    
                    geneJson <- fromJSON(toJSON(content(r)))
                    
                    break()
                    
                } else {
                    
                    headers <- headers(r)
                    
                    remaining <- headers[["x-ratelimit-remaining"]]
                    reset <- headers[["x-ratelimit-reset"]]
                    
                    if (!is.null(remaining) && !is.null(reset) && remaining <= 1) {
                        
                        Sys.sleep(reset + 1)
                        
                    } else {
                        
                        Sys.sleep(try)
                        
                    }
                }
            },
            error = function(e) {
                
                print(e)
                vepJson <- NULL
                
            }
        )
        
    }
    
    if (!is.null(geneJson) && length(geneJson) > 0) {
        
        for (i in 1:length(geneJson)) {
            
            geneId <- geneJson[[i]]$id
            
            if (startsWith(geneId, "ENSG")) {
                
                return(geneId)
                
            }
        }
    }
    
    return(NULL)
    
}

getGeneCoordinates <- function(
    geneIdentifier
) {
    
    ext <- paste0("/overlap/id/", geneIdentifier, "?feature=gene")
    
    for (try in 1:100) {
        
        tryCatch(
            {
                r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
                
                responseCode <- status_code(r)
                
                if (responseCode == 200) {
                    
                    geneJson <- fromJSON(toJSON(content(r)))
                    
                    break()
                    
                } else {
                    
                    headers <- headers(r)
                    
                    remaining <- headers[["x-ratelimit-remaining"]]
                    reset <- headers[["x-ratelimit-reset"]]
                    
                    if (!is.null(remaining) && !is.null(reset) && remaining <= 1) {
                        
                        Sys.sleep(reset + 1)
                        
                    } else {
                        
                        Sys.sleep(try)
                        
                    }
                }
            },
            error = function(e) {
                
                print(e)
                geneJson <- NULL
                
            }
        )
        
    }
    
    if (!is.null(geneJson) && length(geneJson) > 0) {
        
        for (i in 1:length(geneJson)) {
            
            feature <- geneJson[[i]]
            
            if (feature$feature_type == "gene" && feature$gene_id == geneIdentifier) {
                
                chr <- feature$seq_region_name
                start <- feature$start
                end <- feature$end
                
                return(c(chr, start, end))
                
            }
        }
    }
    
    return(c(NA, NA))
    
}

# Files

targetsTable <- "resources/targets/targets_pw_20.01.22"
lmPwTable <- "docs/results/trio_moba/targets_pw_z_placenta_weight_20.01.22.lm.gz"
lmBwTable <- "docs/results/trio_moba/targets_pw_z_weight0_20.01.22.lm.gz"
decodeTable <- "resources/imprinting/decode_clustering.txt"
decodeAlleleTable <- "resources/imprinting/decode_ea.txt"
imprintedGenesTable <- "resources/imprinting/imprinted_genes.txt"
paperTable1 <- "resources/imprinting/table_1_20.01.22.txt"


# Load data

decodeDF <- read.table(
    file = decodeTable,
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>% 
    clean_names()

deacodeAllelesDF <- read.table(
    file = decodeAlleleTable,
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>% 
    clean_names()

imprintedGenesDF <- read.table(
    file = imprintedGenesTable,
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>% 
    clean_names()

targetsDF <- read.table(
    file = targetsTable,
    header = T,
    sep = "\t",
    stringsAsFactors = F
)

lmPwDF <- read.table(
    file = lmPwTable,
    header = T,
    sep = "\t",
    stringsAsFactors = F
)

lmBwDF <- read.table(
    file = lmBwTable,
    header = T,
    sep = "\t",
    stringsAsFactors = F
)

paperTable1DF <- read.table(
    file = paperTable1,
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>% 
    clean_names()


# Add effect allele information to decode table

deacodeAllelesDF <- deacodeAllelesDF %>% 
    group_by(
        variant
    ) %>% 
    arrange(
        p
    ) %>% 
    filter(
        row_number() == 1
    )

decodeDF$ea <- NA
decodeDF$oa <- NA
decodeDF$eaf <- NA

for (i in 1:nrow(decodeDF)) {
    
    j <- which(deacodeAllelesDF$variant == decodeDF$variant[i])
    
    if (deacodeAllelesDF$genome[j] == "Child") {
        
        if (sign(decodeDF$ft_beta[i]) == sign(deacodeAllelesDF$beta[j])) {
            
            swap <- F
            
        } else {
            
            swap <- T
            
        }
        
    } else if (deacodeAllelesDF$genome[j] == "Mother") {
        
        if (sign(decodeDF$mnt_beta[i]) == sign(deacodeAllelesDF$beta[j])) {
            
            swap <- F
            
        } else {
            
            swap <- T
            
        }
    }
    
    if (swap) {
        
        decodeDF$ea[i] <- deacodeAllelesDF$oa[j]
        decodeDF$oa[i] <- deacodeAllelesDF$ea[j]
        decodeDF$eaf[i] <- 1 - deacodeAllelesDF$eaf[j] / 100
        
    } else {
        
        decodeDF$ea[i] <- deacodeAllelesDF$ea[j]
        decodeDF$oa[i] <- deacodeAllelesDF$oa[j]
        decodeDF$eaf[i] <- deacodeAllelesDF$eaf[j] / 100
        
    }
}



# Re-add the rsid for the multi-allelic variant, remove the non-relevant allele

lmPwDF$rsid[lmPwDF$variantId == "5_157808173_A_T"] <- "rs72804545"
lmBwDF$rsid[lmBwDF$variantId == "5_157808173_A_T"] <- "rs72804545"
targetsDF$snp[targetsDF$snp == "5_157808173_A_T"] <- "rs72804545"

lmPwDF <- lmPwDF %>% 
    filter(
        variantId != "5_157808173_A_G"
    )
lmBwDF <- lmBwDF %>% 
    filter(
        variantId != "5_157808173_A_G"
    )


# Get data from the other tables

annotatedTargetsDF <- targetsDF %>% 
    left_join(
        paperTable1DF %>% 
            select(
                snp = rsid,
                maternal_fetal_classification,
                pw_bw_classification
            ),
        by = "snp"
    )


# Imprinted genes coordinates

imprintedGenesDF$chr <- NA
imprintedGenesDF$start <- NA
imprintedGenesDF$end <- NA

for (i in 1:nrow(imprintedGenesDF)) {
    
    geneName <- imprintedGenesDF$gene[i]
    
    geneIdentifier <- getGeneIdentifier(geneName)
    
    if (is.null(geneIdentifier)) {
        
        geneNames <- strsplit(imprintedGenesDF$aliases[i], ",")[[1]]
        
        for (geneName in geneNames) {
            
            geneIdentifier <- getGeneIdentifier(geneName)
            
            if (!is.null(geneIdentifier)) {
                
                break()
                
            }
        }
    }
    
    if (!is.null(geneIdentifier)) {
        
        geneCoordinates <- getGeneCoordinates(geneIdentifier)
        
        imprintedGenesDF$chr[i] <- geneCoordinates[1]
        imprintedGenesDF$start[i] <- geneCoordinates[2]
        imprintedGenesDF$end[i] <- geneCoordinates[3]
        
    }
}

imprintedGenesDF$start <- as.numeric(imprintedGenesDF$start)
imprintedGenesDF$end <- as.numeric(imprintedGenesDF$end)


# Distance to closest imprinted gene

imprintedGeneDistanceDF <- annotatedTargetsDF %>% 
    mutate(
        closest_imprinted_gene = NA,
        closest_imprinted_gene_status = NA,
        closest_imprinted_gene_expressed_allele = NA,
        distance_imprinted_gene = NA
    )

for (i in 1:nrow(targetsDF)) {
    
    bestDistance <- NA
    
    snpChr <- targetsDF$chr[i]
    snpPos <- targetsDF$start[i]
    
    tempDF <- imprintedGenesDF %>% 
        filter(
            !is.na(start) & !is.na(end) & !is.na(chr) & snpChr == chr
        ) %>% 
        mutate(
            distance = case_when(
                snpPos < start ~ start - snpPos,
                snpPos > end ~ snpPos - end,
                TRUE ~ 0
            )
        ) %>% 
        arrange(
            distance
        )
    
    if (nrow(tempDF) > 0) {
        
        imprintedGeneDistanceDF$closest_imprinted_gene[i] <- tempDF$gene[1]
        imprintedGeneDistanceDF$closest_imprinted_gene_status[i] <- tempDF$status[1]
        imprintedGeneDistanceDF$closest_imprinted_gene_expressed_allele[i] <- tempDF$expressed_allele[1]
        imprintedGeneDistanceDF$distance_imprinted_gene[i] <- tempDF$distance[1]
        
    }
    
}


# Add haplotype analysis results

lmResultsDF <- imprintedGeneDistanceDF %>% 
    select(
        rsid = snp, locus = comment, maternal_fetal_classification, pw_bw_classification, closest_imprinted_gene, closest_imprinted_gene_status, closest_imprinted_gene_expressed_allele, distance_imprinted_gene
    ) %>% 
    left_join(
        lmPwDF %>% 
            select(
                contig, position, variantId, rsid, testedAllele, otherAllele, n, nAlt, nH, mendelianError, 
                h.variance_explained, h.intercept.p, h.Bmnt, h.Bmt, h.Bft, h.Bfnt, h.Bmnt.se, h.Bmt.se, h.Bft.se, h.Bfnt.se, h.Bmnt.p, h.Bmt.p, h.Bft.p, h.Bfnt.p,
                cmf.variance_explained, cmf.intercept.p, cmf.h.p, cmf.Bc, cmf.Bm, cmf.Bf, cmf.Bc.se, cmf.Bm.se, cmf.Bf.se, cmf.Bc.p, cmf.Bm.p, cmf.Bf.p, 
                cmf_mt.variance_explained, cmf_mt.intercept.p, cmf_mt.Bc, cmf_mt.Bm, cmf_mt.Bf, cmf_mt.Bmt, cmf_mt.Bc.se, cmf_mt.Bm.se, cmf_mt.Bf.se, cmf_mt.Bmt.se, 
                cmf_mt.Bc.p, cmf_mt.Bm.p, cmf_mt.Bf.p, cmf_mt.Bmt.p
            ) %>% 
            rename(
                h.variance_explained.moba.pw = h.variance_explained,
                h.intercept.p.moba.pw = h.intercept.p, 
                h.Bmnt.moba.pw = h.Bmnt, 
                h.Bmt.moba.pw = h.Bmt, 
                h.Bft.moba.pw = h.Bft, 
                h.Bfnt.moba.pw = h.Bfnt, 
                h.Bmnt.se.moba.pw = h.Bmnt.se, 
                h.Bmt.se.moba.pw = h.Bmt.se, 
                h.Bft.se.moba.pw = h.Bft.se, 
                h.Bfnt.se.moba.pw = h.Bfnt.se, 
                h.Bmnt.p.moba.pw = h.Bmnt.p, 
                h.Bmt.p.moba.pw = h.Bmt.p, 
                h.Bft.p.moba.pw = h.Bft.p, 
                h.Bfnt.p.moba.pw = h.Bfnt.p,
                cmf.variance_explained.moba.pw = cmf.variance_explained, 
                cmf.intercept.p.moba.pw = cmf.intercept.p, 
                cmf.h.p.moba.pw = cmf.h.p, 
                cmf.Bc.moba.pw = cmf.Bc, 
                cmf.Bm.moba.pw = cmf.Bm, 
                cmf.Bf.moba.pw = cmf.Bf, 
                cmf.Bc.se.moba.pw = cmf.Bc.se, 
                cmf.Bm.se.moba.pw = cmf.Bm.se, 
                cmf.Bf.se.moba.pw = cmf.Bf.se, 
                cmf.Bc.p.moba.pw = cmf.Bc.p, 
                cmf.Bm.p.moba.pw = cmf.Bm.p, 
                cmf.Bf.p.moba.pw = cmf.Bf.p, 
                cmf_mt.variance_explained.moba.pw = cmf_mt.variance_explained, 
                cmf_mt.intercept.p.moba.pw = cmf_mt.intercept.p, 
                cmf_mt.Bc.moba.pw = cmf_mt.Bc, 
                cmf_mt.Bm.moba.pw = cmf_mt.Bm, 
                cmf_mt.Bf.moba.pw = cmf_mt.Bf, 
                cmf_mt.Bmt.moba.pw = cmf_mt.Bmt, 
                cmf_mt.Bc.se.moba.pw = cmf_mt.Bc.se, 
                cmf_mt.Bm.se.moba.pw = cmf_mt.Bm.se, 
                cmf_mt.Bf.se.moba.pw = cmf_mt.Bf.se, 
                cmf_mt.Bmt.se.moba.pw = cmf_mt.Bmt.se, 
                cmf_mt.Bc.p.moba.pw = cmf_mt.Bc.p, 
                cmf_mt.Bm.p.moba.pw = cmf_mt.Bm.p, 
                cmf_mt.Bf.p.moba.pw = cmf_mt.Bf.p, 
                cmf_mt.Bmt.p.moba.pw = cmf_mt.Bmt.p
            ),
        by = "rsid"
    ) %>% 
    left_join(
        lmBwDF %>% 
            select(
                rsid, 
                h.variance_explained, h.intercept.p, h.Bmnt, h.Bmt, h.Bft, h.Bfnt, h.Bmnt.se, h.Bmt.se, h.Bft.se, h.Bfnt.se, h.Bmnt.p, h.Bmt.p, h.Bft.p, h.Bfnt.p,
                cmf.variance_explained, cmf.intercept.p, cmf.h.p, cmf.Bc, cmf.Bm, cmf.Bf, cmf.Bc.se, cmf.Bm.se, cmf.Bf.se, cmf.Bc.p, cmf.Bm.p, cmf.Bf.p, 
                cmf_mt.variance_explained, cmf_mt.intercept.p, cmf_mt.Bc, cmf_mt.Bm, cmf_mt.Bf, cmf_mt.Bmt, cmf_mt.Bc.se, cmf_mt.Bm.se, cmf_mt.Bf.se, cmf_mt.Bmt.se, 
                cmf_mt.Bc.p, cmf_mt.Bm.p, cmf_mt.Bf.p, cmf_mt.Bmt.p
            ) %>% 
            rename(
                h.variance_explained.moba.bw = h.variance_explained,
                h.intercept.p.moba.bw = h.intercept.p, 
                h.Bmnt.moba.bw = h.Bmnt, 
                h.Bmt.moba.bw = h.Bmt, 
                h.Bft.moba.bw = h.Bft, 
                h.Bfnt.moba.bw = h.Bfnt, 
                h.Bmnt.se.moba.bw = h.Bmnt.se, 
                h.Bmt.se.moba.bw = h.Bmt.se, 
                h.Bft.se.moba.bw = h.Bft.se, 
                h.Bfnt.se.moba.bw = h.Bfnt.se, 
                h.Bmnt.p.moba.bw = h.Bmnt.p, 
                h.Bmt.p.moba.bw = h.Bmt.p, 
                h.Bft.p.moba.bw = h.Bft.p, 
                h.Bfnt.p.moba.bw = h.Bfnt.p,
                cmf.variance_explained.moba.bw = cmf.variance_explained, 
                cmf.intercept.p.moba.bw = cmf.intercept.p, 
                cmf.h.p.moba.bw = cmf.h.p, 
                cmf.Bc.moba.bw = cmf.Bc, 
                cmf.Bm.moba.bw = cmf.Bm, 
                cmf.Bf.moba.bw = cmf.Bf, 
                cmf.Bc.se.moba.bw = cmf.Bc.se, 
                cmf.Bm.se.moba.bw = cmf.Bm.se, 
                cmf.Bf.se.moba.bw = cmf.Bf.se, 
                cmf.Bc.p.moba.bw = cmf.Bc.p, 
                cmf.Bm.p.moba.bw = cmf.Bm.p, 
                cmf.Bf.p.moba.bw = cmf.Bf.p, 
                cmf_mt.variance_explained.moba.bw = cmf_mt.variance_explained, 
                cmf_mt.intercept.p.moba.bw = cmf_mt.intercept.p, 
                cmf_mt.Bc.moba.bw = cmf_mt.Bc, 
                cmf_mt.Bm.moba.bw = cmf_mt.Bm, 
                cmf_mt.Bf.moba.bw = cmf_mt.Bf, 
                cmf_mt.Bmt.moba.bw = cmf_mt.Bmt, 
                cmf_mt.Bc.se.moba.bw = cmf_mt.Bc.se, 
                cmf_mt.Bm.se.moba.bw = cmf_mt.Bm.se, 
                cmf_mt.Bf.se.moba.bw = cmf_mt.Bf.se, 
                cmf_mt.Bmt.se.moba.bw = cmf_mt.Bmt.se, 
                cmf_mt.Bc.p.moba.bw = cmf_mt.Bc.p, 
                cmf_mt.Bm.p.moba.bw = cmf_mt.Bm.p, 
                cmf_mt.Bf.p.moba.bw = cmf_mt.Bf.p, 
                cmf_mt.Bmt.p.moba.bw = cmf_mt.Bmt.p
            ),
        by = "rsid"
    )


# Map to decode data

mergeDecodeDF <- lmResultsDF %>% 
    mutate(
        ea_decode = NA,
        oa_decode = NA,
        eaf_moba = NA,
        eaf_decode = NA,
        proxy_decode = NA,
        proxy_decode_allele_mapping = NA,
        proxy_decode_r2 = NA,
        group_decode = NA,
        cluster_decode = NA,
        decode_ft_beta = NA,
        decode_ft_se = NA,
        decode_ft_p = NA,
        decode_mt_beta = NA,
        decode_mt_se = NA,
        decode_mt_p = NA,
        decode_mnt_beta = NA,
        decode_mnt_se = NA,
        decode_mnt_p = NA
    )

proxies <- list(
    rs72804545 = "rs78190170"
)


for (i in 1:nrow(mergeDecodeDF)) {
    
    variant <- mergeDecodeDF$rsid[i]
    
    if (variant %in% names(proxies)) {
        
        variant <- proxies[[variant]]
        
    }
    
    eafText <- mergeDecodeDF$nAlt[i]
    
    eafMother <- strsplit(strsplit(eafText, "mother\\(")[[1]][2], "\\)")[[1]][1]
    eafFather <- strsplit(strsplit(eafText, "father\\(")[[1]][2], "\\)")[[1]][1]
    
    eafMother <- strsplit(eafMother, ",")[[1]]
    eafFather <- strsplit(eafFather, ",")[[1]]
    
    n0Mother <- as.numeric(strsplit(eafMother[1], ":")[[1]][2])
    n0Father <- as.numeric(strsplit(eafFather[1], ":")[[1]][2])
    n1Mother <- as.numeric(strsplit(eafMother[2], ":")[[1]][2])
    n1Father <- as.numeric(strsplit(eafFather[2], ":")[[1]][2])
    
    if (length(eafMother) == 3) {
        
        n2Mother <- as.numeric(strsplit(eafMother[3], ":")[[1]][2])
        
    } else {
        
        n2Mother <- 0
        
    }
    
    if (length(eafFather) == 3) {
        
        n2Father <- as.numeric(strsplit(eafFather[3], ":")[[1]][2])
        
    } else {
        
        n2Father <- 0
        
    }
    
    n0 <- n0Mother + n0Father
    n1 <- n1Mother + n1Father
    n2 <- n2Mother + n2Father
    
    
    eafMoba <- (2* n2 + n1) / (2*(n0 + n1 + n2))
    
    mergeDecodeDF$eaf_moba[i] <- eafMoba
    
    j <- -1
    
    if (variant %in% decodeDF$variant) {
        
        j <- which(decodeDF$variant == variant)
        
        if (decodeDF$ea[j] == mergeDecodeDF$otherAllele[i] && decodeDF$oa[j] == mergeDecodeDF$testedAllele[i]) {
            
            decodeEa <- decodeDF$oa[j]
            decodeOa <- decodeDF$ea[j]
            decodeEaf <- 1 - decodeDF$eaf[j]
            betaCoefficient <- -1
            
        } else if (decodeDF$ea[j] == mergeDecodeDF$testedAllele[i] && decodeDF$oa[j] == mergeDecodeDF$otherAllele[i]) {
            
            decodeEa <- decodeDF$ea[j]
            decodeOa <- decodeDF$oa[j]
            decodeEaf <- decodeDF$eaf[j]
            swap <- F
            betaCoefficient <- 1
            
        } else {
            
            stop(glue("Allele mismatch for {variant}: decode {decodeDF$ea[j]},{decodeDF$oa[j]}; {mergeDecodeDF$testedAllele[i]},{mergeDecodeDF$otherAllele[i]}."))
            
        }
        
        mergeDecodeDF$ea_decode[i] <- decodeEa
        mergeDecodeDF$oa_decode[i] <- decodeOa
        mergeDecodeDF$eaf_decode[i] <- decodeEaf
        mergeDecodeDF$group_decode[i] <- decodeDF$group[j]
        mergeDecodeDF$cluster_decode[i] <- decodeDF$cluster[j]
        mergeDecodeDF$decode_ft_beta[i] <- betaCoefficient * decodeDF$ft_beta[j]
        mergeDecodeDF$decode_ft_se[i] <- decodeDF$ft_se[j]
        mergeDecodeDF$decode_ft_p[i] <- decodeDF$ft_p[j]
        mergeDecodeDF$decode_mt_beta[i] <- betaCoefficient * decodeDF$mt_beta[j]
        mergeDecodeDF$decode_mt_se[i] <- decodeDF$mt_se[j]
        mergeDecodeDF$decode_mt_p[i] <- decodeDF$mt_p[j]
        mergeDecodeDF$decode_mnt_beta[i] <- betaCoefficient * decodeDF$mnt_beta[j]
        mergeDecodeDF$decode_mnt_se[i] <- decodeDF$mnt_se[j]
        mergeDecodeDF$decode_mnt_p[i] <- decodeDF$mnt_p[j]
        
    } else {
        
        proxy <- LDproxy(
            snp = variant, 
            pop = "CEU", 
            token = "972f33fe5966"
        )
        
        proxy <- proxy %>% 
            clean_names() %>% 
            mutate(
                rs_number = as.character(rs_number),
                coord = as.character(coord),
                alleles = as.character(alleles),
                correlated_alleles = as.character(correlated_alleles)
            ) %>% 
            filter(
                rs_number %in% decodeDF$variant
            ) %>% 
            arrange(
                desc(r2)
            )
        
        if (nrow(proxy) > 0) {
            
            j <- which(decodeDF$variant == proxy$rs_number[1])
        
            correlatedAlleles <- proxy$correlated_alleles[1]
            
            alleleMatch1 <- glue("{mergeDecodeDF$testedAllele[i]}={decodeDF$ea[j]},{mergeDecodeDF$otherAllele[i]}={decodeDF$oa[j]}")
            alleleMatch2 <- glue("{mergeDecodeDF$otherAllele[i]}={decodeDF$oa[j]},{mergeDecodeDF$testedAllele[i]}={decodeDF$ea[j]}")
            alleleMatch3 <- glue("{decodeDF$ea[j]}={mergeDecodeDF$testedAllele[i]},{decodeDF$oa[j]}={mergeDecodeDF$otherAllele[i]}")
            alleleMatch4 <- glue("{decodeDF$oa[j]}={mergeDecodeDF$otherAllele[i]},{decodeDF$ea[j]}={mergeDecodeDF$testedAllele[i]}")
            alleleSwap1 <- glue("{mergeDecodeDF$testedAllele[i]}={decodeDF$oa[j]},{mergeDecodeDF$otherAllele[i]}={decodeDF$ea[j]}")
            alleleSwap2 <- glue("{mergeDecodeDF$otherAllele[i]}={decodeDF$ea[j]},{mergeDecodeDF$testedAllele[i]}={decodeDF$oa[j]}")
            alleleSwap3 <- glue("{decodeDF$oa[j]}={mergeDecodeDF$testedAllele[i]},{decodeDF$ea[j]}={mergeDecodeDF$otherAllele[i]}")
            alleleSwap4 <- glue("{decodeDF$ea[j]}={mergeDecodeDF$otherAllele[i]},{decodeDF$oa[j]}={mergeDecodeDF$testedAllele[i]}")
            
            decodeEaf <- decodeDF$eaf[j]
            mobaEaf <- mergeDecodeDF$eaf_moba[i]
            
            if (is.na(correlatedAlleles)) {
                
                if (decodeEaf <= 0.5 && mobaEaf > 0.5 || decodeEaf > 0.5 && mobaEaf <= 0.5) {
                    
                    proxySwap <- T
                    
                } else {
                    
                    proxySwap <- F
                    
                }
                
            } else if (correlatedAlleles == alleleMatch1 || correlatedAlleles == alleleMatch2 || correlatedAlleles == alleleMatch3 || correlatedAlleles == alleleMatch4) {
                
                proxySwap <- F
                
            } else if (correlatedAlleles == alleleSwap1 || correlatedAlleles == alleleSwap2 || correlatedAlleles == alleleSwap3 || correlatedAlleles == alleleSwap4) {
                
                proxySwap <- T
                
            } else {
                
                stop(glue("Allele combination not supported: '{correlatedAlleles}', decode '{decodeDF$ea[j]},{decodeDF$oa[j]}', moba '{mergeDecodeDF$testedAllele[i]},{mergeDecodeDF$otherAllele[i]}'"))
                
            }
            
            if (proxySwap) {
                
                betaSign <- -1
                decodeEa <- decodeDF$oa[j]
                decodeOa <- decodeDF$ea[j]
                decodeEaf <- 1 - decodeEaf
                
            } else {
                
                betaSign <- 1
                decodeEa <- decodeDF$ea[j]
                decodeOa <- decodeDF$oa[j]
                
            }
            
            mergeDecodeDF$ea_decode[i] <- decodeEa
            mergeDecodeDF$oa_decode[i] <- decodeOa
            mergeDecodeDF$eaf_decode[i] <- decodeEaf
            mergeDecodeDF$proxy_decode[i] <- proxy$rs_number[1]
            mergeDecodeDF$proxy_decode_allele_mapping[i] <- proxy$correlated_alleles[1]
            mergeDecodeDF$proxy_decode_r2[i] <- proxy$r2[1]
            mergeDecodeDF$group_decode[i] <- decodeDF$group[j]
            mergeDecodeDF$cluster_decode[i] <- decodeDF$cluster[j]
            mergeDecodeDF$decode_ft_beta[i] <- betaSign * decodeDF$ft_beta[j]
            mergeDecodeDF$decode_ft_se[i] <- decodeDF$ft_se[j]
            mergeDecodeDF$decode_ft_p[i] <- decodeDF$ft_p[j]
            mergeDecodeDF$decode_mt_beta[i] <- betaSign * decodeDF$mt_beta[j]
            mergeDecodeDF$decode_mt_se[i] <- decodeDF$mt_se[j]
            mergeDecodeDF$decode_mt_p[i] <- decodeDF$mt_p[j]
            mergeDecodeDF$decode_mnt_beta[i] <- betaSign * decodeDF$mnt_beta[j]
            mergeDecodeDF$decode_mnt_se[i] <- decodeDF$mnt_se[j]
            mergeDecodeDF$decode_mnt_p[i] <- decodeDF$mnt_p[j]
            
        }
    }
}

write.table(
    x = mergeDecodeDF,
    file = "docs/figures/trio/merge_decode_27.01.22.txt",
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t"
)



