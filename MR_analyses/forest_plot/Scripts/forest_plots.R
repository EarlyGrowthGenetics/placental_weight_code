library(here)
library(glue)
library(dplyr)
library(ggplot2)
library(grid)

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

df$y <- as.numeric(df$person)
df$y <- ifelse(df$y == min(df$y), df$y + 0.5, df$y)
df$y <- ifelse(df$y == max(df$y), df$y - 0.5, df$y)


#### Change column names for legend

df$Analysis[df$Analysis == "Meta"] <- "Meta-analysis"
df$Analysis[df$Analysis == "WLM"] <- "WLM-adjusted"


###### Plot

forest_plot <- ggplot() + 
  geom_vline(
    xintercept = 0
  ) +
  geom_segment(
    data = df,
    mapping = aes(
      x = Beta - qnorm(0.975) * Se,
      xend = Beta + qnorm(0.975) * Se,
      y = y,
      yend = y,
      color = Genome
    ),
    size = 0.8
  ) +
  geom_point(
    data = df,
    mapping = aes(
      x = Beta,
      y = y,
      color = Genome,
      shape = Analysis
    ),
    size = 4
  ) + 
  facet_grid(
    exposure ~ .,
    scales = "free",
    space = "free",
    labeller = label_wrap_gen(
      width = 2,
      multi_line = T
      )
    ) +
  scale_x_continuous(
    name = "Effect Size [95% CI]"
  ) +
  scale_y_continuous(
    breaks = c(1.75, 3.25),
    labels = c("Fetal", "Maternal"),
    limits = c(1, 4)
  ) +
  scale_color_manual(
    values = c("#009E73","#D55E00"),
    guide = "none"
    ) + 
  theme_bw(
      base_size = 24
    ) + 
  scale_shape(
    guide = guide_legend(
      override.aes = list(size = 5)
    )
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y.left = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.position = c(0.13, 0.96),
    legend.title = element_blank(),
    legend.background = element_rect(
      fill = "white",
      color = "grey20",
      size = 0.5
    ),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 18, face = "bold"),
    legend.margin = margin(t = -3, r = 2, b = 2, l = 2, unit = "mm")
  ) 

# save the plot

figure_path <- file.path(here(), 'MR_analyses/forest_plot/Fig6_Mendelian_Randomization.png')
png(figure_path, width = 900, height = 900)
grid.draw(forest_plot)
device <- dev.off()

