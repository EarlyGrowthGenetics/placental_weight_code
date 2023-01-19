
##
#
# This script makes a forest plot of the association results stratified by parity.
#
##


# Packages
library(here)
library(glue)
library(tidyr)
library(dplyr)
library(janitor)


# Load variants
variants_table <- read.table(
  file = file.path(here(), "parity/targets_pw_20.01.22"),
  sep = "\t",
  header = T,
  stringsAsFactors = F
)


# Load trios
sample_list <- read.table(
  file = "/mnt/archive/marc/back-up_Harvest/trio/trio",
  sep = "\t",
  header = T,
  stringsAsFactors = F
)


# Load genotypes
genotypes <- list()

for (chr in unique(variants_table$chr)) {
  
  genotypes_table <- read.table(
    file = glue("/mnt/work/marc/egg/placental_weight/data/genotypes_top_snps/", chr, ".raw"),
    sep = "\t",
    header = T,
    stringsAsFactors = F
  ) %>% 
    clean_names() %>% 
    select(
      iid, starts_with("rs")
    ) %>% 
    pivot_longer(
      cols = starts_with("rs"),
      names_to = "variant"
    )
  
  genotypes[[length(genotypes) + 1]] <- genotypes_table
  
}

genotypes <- do.call("rbind", genotypes)


# Load PCs

pcs_table <- read.table(
  file = "/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/pca/mobagen-total/mobagen-total-proj-pc",
  header = T,
  sep = "\t",
  stringsAsFactors = F
) %>% 
  clean_names()


# Load phenotypes
pregnancy_table <- read.table(
  file = "/mnt/archive/moba/pheno/v10/V10_1.1.1-221121/pregnancy.gz",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% 
  select(
    child_id, sex, pregnancy_duration, plural_birth, rank_siblings, n_previous_deliveries
  )

delivery_table <- read.table(
  file = "/mnt/archive/moba/pheno/v10/V10_1.1.1-221121/delivery.gz",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% 
  select(
    child_id, child_SentrixID, mother_SentrixID, father_SentrixID, placenta_weight, birth_defect
  )


# Exclusion criteria

pheno_table <- pregnancy_table %>% 
  full_join(
    delivery_table,
    by = "child_id"
  ) %>% 
  filter(
    !is.na(pregnancy_duration) & !is.na(placenta_weight) & pregnancy_duration > 37 * 7 & pregnancy_duration < 43 * 7 & plural_birth == "Single birth" & is.na(birth_defect) & placenta_weight >= 200 & placenta_weight <= 1500
  ) %>% 
  mutate(
    parity = ifelse(n_previous_deliveries == 0, 0, 2),
    parity = ifelse(n_previous_deliveries == 1, 1, parity),
  )


# Make individual level DFs

child_table <- pheno_table %>% 
  filter(
    child_SentrixID %in% sample_list$child_SentrixID
  ) %>% 
  select(
    iid = child_SentrixID, sex, pregnancy_duration, parity, placenta_weight
  ) %>% 
  left_join(
    pcs_table %>% 
      select(
        iid, starts_with("pc")
      ),
    by = "iid"
  ) %>% 
  mutate(
    z_placenta_weight = (placenta_weight - mean(placenta_weight)) / sd(placenta_weight),
    z_parity = (parity - mean(parity)) / sd(parity)
  ) %>% 
  filter(
    abs(z_placenta_weight) <= 5
  )


# Run linear model with parity

parity_lm_coefficients <- list()
pw_lm_coefficients <- list()
stratified_pw_lm_coefficients <- list()

for (variant_id in variants_table$snp) {
  
  snp_table <- genotypes %>% 
    filter(
      startsWith(variant, variant_id)
    ) %>% 
    inner_join(
      child_table,
      by = "iid"
    )
  
  if (nrow(snp_table) > 0) {
    
    # Parity
    
    lm_results <- lm(
      formula = "z_parity ~ value + sex + pregnancy_duration + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10",
      data = snp_table
    )
    
    lm_results_summary <- summary(lm_results)
    
    lm_results_summary_coefficients <- as.data.frame(lm_results_summary$coefficients)
    
    lm_results_summary_coefficients$snp <- variant_id
    
    parity_lm_coefficients[[length(parity_lm_coefficients) + 1]] <- lm_results_summary_coefficients
    
    
    # PW
    
    lm_results <- lm(
      formula = "z_placenta_weight ~ value + sex + pregnancy_duration + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10",
      data = snp_table
    )
    
    lm_results_summary <- summary(lm_results)
    
    lm_results_summary_coefficients <- as.data.frame(lm_results_summary$coefficients)
    
    lm_results_summary_coefficients$snp <- variant_id
    lm_results_summary_coefficients$model <- "pw"
    
    pw_lm_coefficients[[length(pw_lm_coefficients) + 1]] <- lm_results_summary_coefficients
    
    
    # PW + parity
    
    lm_results <- lm(
      formula = "z_placenta_weight ~ value + sex + z_parity + pregnancy_duration + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10",
      data = snp_table
    )
    
    lm_results_summary <- summary(lm_results)
    
    lm_results_summary_coefficients <- as.data.frame(lm_results_summary$coefficients)
    
    lm_results_summary_coefficients$snp <- variant_id
    lm_results_summary_coefficients$model <- "pw + parity"
    
    pw_lm_coefficients[[length(pw_lm_coefficients) + 1]] <- lm_results_summary_coefficients
    
    
    # Stratify
    
    for (parity_level in 0:2) {
      
      temp_table <- snp_table %>% 
        filter(
          parity == parity_level
        )
      
      lm_results <- lm(
        formula = "z_placenta_weight ~ value + sex + pregnancy_duration + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10",
        data = temp_table
      )
      
      lm_results_summary <- summary(lm_results)
      
      lm_results_summary_coefficients <- as.data.frame(lm_results_summary$coefficients)
      
      lm_results_summary_coefficients$snp <- variant_id
      lm_results_summary_coefficients$parity_level <- parity_level
      
      stratified_pw_lm_coefficients[[length(stratified_pw_lm_coefficients) + 1]] <- lm_results_summary_coefficients
      
    }
  }
}

parity_lm_coefficients <- do.call("rbind", parity_lm_coefficients)
pw_lm_coefficients <- do.call("rbind", pw_lm_coefficients)
stratified_pw_lm_coefficients <- do.call("rbind", stratified_pw_lm_coefficients)

parity_lm_coefficients$variable <- row.names(parity_lm_coefficients)
pw_lm_coefficients$variable <- row.names(pw_lm_coefficients)
stratified_pw_lm_coefficients$variable <- row.names(stratified_pw_lm_coefficients)


# Export results

write.table(
  parity_lm_coefficients,
  file = gzfile(file.path(here(), "parity/lm/parity_lm_coefficients.gz")),
  col.names = T,
  row.names = F,
  sep = "\t"
)

write.table(
  pw_lm_coefficients,
  file = gzfile(file.path(here(), "parity/lm/pw_lm_coefficients.gz")),
  col.names = T,
  row.names = F,
  sep = "\t"
)

write.table(
  stratified_pw_lm_coefficients,
  file = gzfile(file.path(here(), "parity/lm/stratified_pw_lm_coefficients.gz")),
  col.names = T,
  row.names = F,
  sep = "\t"
)





