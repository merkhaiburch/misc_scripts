# Subset by trait, remove all 0 p values
gwa_results <- read.delim("~/Desktop/glm_maxdiv0.001.txt")

# Subset dataframe
dpoll <- gwa_results[which(gwa_results$Trait == "dpoll"),]
EarDia <- gwa_results[which(gwa_results$Trait == "EarDia"),]
EarHT <- gwa_results[which(gwa_results$Trait == "EarHT"),]

# Manhattan Plot
source('~/git_projects/haplotype_v_snp_gwas/src/R/random_scripts/manhattan_plot.R')
library(ggplot2)
a <- manhattan_plot(dpoll, sig_threshold = 3, ylim = 5, title = "dpol")
c <- manhattan_plot(EarDia, sig_threshold = 3, ylim = 5, title = "EarDia")
e <- manhattan_plot(EarHT, sig_threshold = 3, ylim = 5, title = "EarHt")

# Setup data
dpoll$observed <- -log10(dpoll$p)
dpoll$expected <- -log10(ppoints(length(dpoll$p)))

EarDia$observed <- -log10(EarDia$p)
EarDia$expected <- -log10(ppoints(length(EarDia$p)))

EarHT$observed <- -log10(EarHT$p)
EarHT$expected <- -log10(ppoints(length(EarHT$p)))

# QQ Plot
b <- ggplot(dpoll, aes(x=expected, y=sort(observed,decreasing = T))) +
  geom_point(size = 1)  +
  xlim(0 , 5) + 
  geom_abline(intercept = 0) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 17)) +
  scale_size_continuous(range = c(0,3)) +
  labs(x = "Observed", y = "Expected") +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 9, vjust = 0.5),
    axis.text.y = element_text(size = 9, vjust = 0.5))

d <- ggplot(EarDia, aes(x=expected, y=sort(observed,decreasing = T))) +
  geom_point(size = 1)  +
  xlim(0 , 5) + 
  geom_abline(intercept = 0) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 17)) +
  scale_size_continuous(range = c(0,3)) +
  labs(x = "Observed", y = "Expected") +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 9, vjust = 0.5),
    axis.text.y = element_text(size = 9, vjust = 0.5))

f <- ggplot(EarHT, aes(x=expected, y=sort(observed,decreasing = T))) +
  geom_point(size = 1)  +
  xlim(0 , 5) + 
  geom_abline(intercept = 0) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 17)) +
  scale_size_continuous(range = c(0,3)) +
  labs(x = "Observed", y = "Expected") +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 9, vjust = 0.5),
    axis.text.y = element_text(size = 9, vjust = 0.5))

# Use patchwork
library(patchwork)
(a+b)/(c+d)/(e+f)

a/c/e

