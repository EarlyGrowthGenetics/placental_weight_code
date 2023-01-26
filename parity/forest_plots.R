
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
library(ggplot2)
library(ggrepel)
library(grid)

# Functions
getExpectedPvalue <- function(observedPValues) {
  
  n <- length(observedPValues)
  
  return(
    ppoints(n = n)
  )
  
}

getCiUp <- function(observedPValues, confidence) {
  
  n <- length(observedPValues)
  
  return(
    qbeta(p = (1-confidence)/2, shape1 = 1:n, shape2 = n - 1:n + 1)
  )
}

getCiDown <- function(observedPValues, confidence) {
  
  n <- length(observedPValues)
  
  return(
    qbeta(p = 1-(1-confidence)/2, shape1 = 1:n, shape2 = n - 1:n + 1)
  )
}

getQqPlot <- function(
    pValues,
    variantIds,
    categories1 = "",
    categories2 = "",
    confidence = 0.95
) {
  
  qqDF <- data.frame(
    pValue = pValues,
    variantId = variantIds,
    category1 = categories1,
    category2 = categories2,
    stringsAsFactors = F
  ) %>%
    group_by(
      category1, category2
    ) %>%
    arrange(
      pValue
    ) %>%
    mutate(
      expectedPValue = getExpectedPvalue(
        observedPValues = pValue
      ),
      ciUp = getCiUp(
        observedPValues = pValue, 
        confidence = confidence
      ),
      ciDown = getCiDown(
        observedPValues = pValue, 
        confidence = confidence
      ),
      logP = -log10(pValue),
      expectedLogP = -log10(expectedPValue),
      ciUpLog = -log10(ciDown),
      ciDownLog = -log10(ciUp)
    ) %>%
    filter(
      !is.na(logP)
    ) %>%
    mutate(
      label = ifelse(row_number() <= 3, variantId, "")
    ) %>% 
    ungroup()
  
  qqPlot <- ggplot() +
    theme_bw(
      base_size = 24
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      col = "black",
      linetype = "dashed"
    ) +
    geom_ribbon(
      data = qqDF,
      mapping = aes(
        x = expectedLogP,
        ymin = ciDownLog,
        ymax = ciUpLog
      ),
      fill = "black",
      alpha = 0.2
    ) +
    geom_point(
      data = qqDF,
      mapping = aes(
        x = expectedLogP,
        y = logP
      )
    ) +
    geom_label_repel(
      data = qqDF,
      mapping = aes(
        x = expectedLogP,
        y = logP,
        label = label
      ),
      point.padding = 0.5,
      nudge_y = 0.5 + qqDF$logP / 10
    ) +
    scale_x_continuous(
      name = "Expected p-value [-log10]"
    ) + 
    scale_y_continuous(
      name = "Observed p-value [-log10]"
    ) + 
    theme(
      legend.position = "none"
    )
  
  if (length(categories1) > 1 && length(categories2) == 1) {
    
    qqPlot <- qqPlot +
      facet_grid(
        . ~ category1
      )
    
  } else if (length(categories1) == 1 && length(categories2) > 1) {
    
    qqPlot <- qqPlot +
      facet_grid(
        category2 ~ .
      )
    
  } else if (length(categories1) > 1 && length(categories2) > 1) {
    
    qqPlot <- qqPlot +
      facet_grid(
        category2 ~ category1
      )
    
  }
  
  return(qqPlot)
  
}



# Load the results of the analysis
variants_table <- read.table(
  file = file.path(here(), "parity/extract_supp_table_7_23.01.18"),
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% 
  clean_names()

parity_lm_coefficients <- read.table(
  file = file.path(here(), "parity/lm/parity_lm_coefficients.gz"),
  header = T,
  sep = "\t",
  stringsAsFactors = F
) %>% 
  clean_names() %>% 
  mutate(
    individual = case_when(
      individual == "child" ~ "Child",pain
      individual == "mother" ~ "Mother",
      individual == "father" ~ "Father"
    )
  )

pw_lm_coefficients <- read.table(
  file = file.path(here(), "parity/lm/pw_lm_coefficients.gz"),
  header = T,
  sep = "\t",
  stringsAsFactors = F
) %>% 
  clean_names() %>% 
  mutate(
    individual = case_when(
      individual == "child" ~ "Child",
      individual == "mother" ~ "Mother",
      individual == "father" ~ "Father"
    )
  )

stratified_pw_lm_coefficients <- read.table(
  file = file.path(here(), "parity/lm/stratified_pw_lm_coefficients.gz"),
  header = T,
  sep = "\t",
  stringsAsFactors = F
) %>% 
  clean_names() %>% 
  mutate(
    individual = case_when(
      individual == "child" ~ "Child",
      individual == "mother" ~ "Mother",
      individual == "father" ~ "Father"
    )
  )


# QQ plot of parity

parity_plot_data <- parity_lm_coefficients %>% 
  filter(
    variable == "value"
  ) %>% 
  select(
    snp,
    beta = estimate,
    se = std_error,
    p = pr_t,
    individual
  ) %>% 
  mutate(
    model = "parity"
  )

qqPlot <- getQqPlot(
  pValues = parity_plot_data$p,
  variantIds = parity_plot_data$snp,
  categories1 = parity_plot_data$individual
)

png(
  filename = file.path(here(), "parity/plots/parity_qq.png"),
  height = 600,
  width = 900
)
grid.draw(qqPlot)
device <- dev.off()


# QQ plot for the interaction term

parity_interaction_plot_data <- pw_lm_coefficients %>% 
  filter(
    variable == "value:z_parity"
  ) %>% 
  select(
    snp,
    beta = estimate,
    se = std_error,
    p = pr_t,
    individual
  )

qqPlotInteraction <- getQqPlot(
  pValues = parity_interaction_plot_data$p,
  variantIds = parity_interaction_plot_data$snp,
  categories1 = parity_interaction_plot_data$individual
)

png(
  filename = file.path(here(), "parity/plots/parity_interaction_qq.png"),
  height = 600,
  width = 900
)
grid.draw(qqPlotInteraction)
device <- dev.off()


# Forest plot for covariates

covariates_plot_data <- pw_lm_coefficients %>% 
  filter(
    variable == "value"
  ) %>% 
  select(
    snp,
    beta = estimate,
    se = std_error,
    p = pr_t,
    covariate = model,
    individual
  )

snp_order_beta <- covariates_plot_data %>% 
  filter(
    covariate == "sex + ga" & individual == "Child"
  ) %>% 
  arrange(abs(beta), p) %>% 
  pull(snp)

snp_order_beta <- rev(snp_order_beta)

snp_order_p <- covariates_plot_data %>% 
  filter(
    covariate == "sex + ga" & individual == "Child"
  ) %>% 
  arrange(p, abs(beta)) %>% 
  pull(snp)

snp_order_p <- rev(snp_order_p)

plot_data <- covariates_plot_data %>% 
  filter(
    covariate %in% c("none", "sex", "sex + ga", "sex + ga + parity")
  ) %>% 
  mutate(
    snp_factor = factor(snp, levels = snp_order_beta),
    covariate_factor = factor(covariate, levels = c("none", "sex", "sex + ga", "sex + ga + parity")),
    y = as.numeric(snp_factor) - 0.2 * (as.numeric(covariate_factor) - 2.5)
  ) %>% 
  arrange(y)

covariates_plot <- ggplot() + 
  theme_bw(
    base_size = 24
  ) +
  geom_vline(
    xintercept = 0
  ) +
  geom_segment(
    data = plot_data,
    mapping = aes(
      x = abs(beta) - qnorm(0.975) * se,
      xend = abs(beta) + qnorm(0.975) * se,
      y = y,
      yend = y,
      col = covariate_factor
    )
  ) +
  geom_point(
    data = plot_data,
    mapping = aes(
      x = abs(beta),
      y = y,
      col = covariate_factor
    )
  ) +
  scale_x_continuous(
    name = "Beta [95% CI]"
  ) +
  scale_y_continuous(
    breaks = 1:length(snp_order_beta),
    labels = snp_order_beta,
    expand = expansion(mult = 0.02)
  ) +
  scale_color_manual(
    name = "Covariates",
    values = c("black", "green3", "blue3", "orange3")
  ) +
  theme(
    legend.position = "top",
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  facet_grid(
    . ~ individual
  )

png(
  filename = file.path(here(), "parity/plots/covariates_effect_size.png"),
  height = 1200,
  width = 900
)
grid.draw(covariates_plot)
device <- dev.off()

plot_data_p <- plot_data %>% 
  mutate(
    snp_factor = factor(snp, levels = snp_order_p),
    y = as.numeric(snp_factor) - 0.2 * (as.numeric(covariate_factor) - 2.5)
  )

covariates_plot_p <- ggplot() +
  theme_bw(
    base_size = 24
  ) +
  geom_vline(
    xintercept = 0
  )  +
  geom_point(
    data = plot_data_p,
    mapping = aes(
      x = -log10(p),
      y = y,
      col = covariate_factor
    )
  ) +
  scale_x_continuous(
    name = "p-value [-log10]"
  ) +
  scale_y_continuous(
    breaks = 1:length(snp_order_p),
    labels = snp_order_p,
    expand = expansion(mult = 0.02)
  ) +
  scale_color_manual(
    name = "Covariates",
    values = c("black", "green3", "blue3", "orange3")
  ) +
  theme(
    legend.position = "top",
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  facet_grid(
    . ~ individual
  )

png(
  filename = file.path(here(), "parity/plots/covariates_p.png"),
  height = 1200,
  width = 900
)
grid.draw(covariates_plot_p)
device <- dev.off()


plot_data_scatter <- covariates_plot_data %>% 
  filter(
    covariate %in% c("sex + ga", "sex + ga + parity")
  ) %>% 
  mutate(
    covariate = ifelse(covariate == "sex + ga", "sex_ga", "sex_ga_parity")
  ) %>% 
  select(
    snp, covariate, individual, p
  ) %>% 
  pivot_wider(
    names_from = covariate,
    values_from = p
  )

covariates_plot_p_scatter <- ggplot() +
  theme_bw(
    base_size = 24
  ) +
  geom_abline(
    linetype = "dotted"
  ) +
  geom_point(
    data = plot_data_scatter,
    mapping = aes(
      x = -log10(sex_ga),
      y = -log10(sex_ga_parity)
    ),
    col = "darkblue",
    alpha = 0.8
  ) +
  scale_x_continuous(
    name = "sex + ga\np-value [-log10]",
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name = "sex + ga + parity\np-value [-log10]",
    expand = expansion(mult = c(0, 0.02))
  ) +
  facet_grid(
    . ~ individual
  )

covariates_plot_p_scatter

png(
  filename = file.path(here(), "parity/plots/covariates_scatter_p.png"),
  height = 600,
  width = 900
)
grid.draw(covariates_plot_p_scatter)
device <- dev.off()

# Forest plot stratified by parity

all_plot_data <- pw_lm_coefficients %>% 
  filter(
    variable == "value" & model == "sex + ga"
  ) %>% 
  select(
    snp,
    beta = estimate,
    se = std_error,
    p = pr_t,
    individual
  ) %>% 
  mutate(
    parity = "all"
  )

stratified_plot_data <- stratified_pw_lm_coefficients %>% 
  filter(
    variable == "value"
  ) %>% 
  select(
    snp,
    beta = estimate,
    se = std_error,
    p = pr_t,
    parity = parity_level,
    individual
  )

plot_data <- rbind(all_plot_data, stratified_plot_data) %>% 
  mutate(
    snp_factor = factor(snp, levels = snp_order_p),
    parity = ifelse(parity == 2, "≥2", parity),
    parity_factor = factor(parity, levels = c("0", "1", "≥2", "all")),
    y = as.numeric(snp_factor) - 0.15 * (as.numeric(parity_factor) - 2.5)
  ) %>% 
  arrange(y)

parity_plot <- ggplot() + 
  theme_bw(
    base_size = 24
  ) +
  geom_vline(
    xintercept = 0
  ) +
  geom_segment(
    data = plot_data,
    mapping = aes(
      x = abs(beta) - qnorm(0.975) * se,
      xend = abs(beta) + qnorm(0.975) * se,
      y = y,
      yend = y,
      col = parity_factor
    )
  ) +
  geom_point(
    data = plot_data,
    mapping = aes(
      x = abs(beta),
      y = y,
      col = parity_factor
    )
  ) +
  scale_x_continuous(
    name = "Beta [95% CI]"
  ) +
  scale_y_continuous(
    breaks = 1:length(snp_order_p),
    labels = snp_order_p,
    expand = expansion(mult = 0.02)
  ) +
  scale_color_manual(
    name = "Parity",
    values = c("darkred", "darkblue", "darkgreen", "black")
  ) +
  theme(
    legend.position = "top",
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  facet_grid(
    . ~ individual
  )

png(
  filename = file.path(here(), "parity/plots/pw_parity.png"),
  height = 900,
  width = 600
)
grid.draw(parity_plot)
device <- dev.off()

