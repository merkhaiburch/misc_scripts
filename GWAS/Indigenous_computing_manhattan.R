# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-04-27
# Updated... 2023-04-27
#
# Description:
# Make manhattan plot of GLM and MLM
# ------------------------------------------------------------------------------

# Load in packages
library(ggplot2)
library(dplyr)

# Load in data
glm <- read.delim("~/Downloads/glm_ear_height.txt")
mlm <- read.delim("~/Downloads/mlm_ear_height.txt")

# all traits: 
unique(glm$Trait)

# Subset data to just ear height
glm <- glm %>% filter(Trait == "EarDia")
mlm <- mlm %>% filter(Trait == "EarDia")

# Format x axis to plot by chromosome cumulative--------------------------------
glm_data_cumulative <- glm %>% 
  group_by(Chr) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chr, bp_add)
glm_gwas_result <- glm %>% 
  inner_join(glm_data_cumulative, by = "Chr") %>% 
  mutate(bp_cumulative = Pos + bp_add)

mlm_data_cumulative <- mlm %>% 
  group_by(Chr) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chr, bp_add)
mlm_gwas_result <- mlm %>% 
  inner_join(mlm_data_cumulative, by = "Chr") %>% 
  mutate(bp_cumulative = Pos + bp_add)


# Plot Manhattan plots ---------------------------------------------------------

# General linear model without k
glm_plot <- ggplot(glm_gwas_result, aes(x = bp_cumulative, y = -log10(p), color = as.factor(Chr), size = -log10(p))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 2.5, color = "grey40", linetype = "solid") +
  labs(x = NULL, y = "-log10(p)") + 
  theme_minimal() +
  theme(legend.position = "none")+
  ylim(0,7)
glm_plot

# Mixed model with K
mlm_plot <- ggplot(mlm_gwas_result, aes(x = bp_cumulative, y = -log10(p), 
                                        color = as.factor(Chr), size = -log10(p))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 2.5, color = "grey40", linetype = "solid") +
  labs(x = NULL, y = "-log10(p)") + 
  theme_minimal() +
  theme(legend.position = "none")+
  ylim(0,7)
mlm_plot


# Export to file ---------------------------------------------------------------

ggsave(plot = glm_plot, filename = "~/Downloads/glm_results.png", width = 10, height = 5, units = "in")
ggsave(plot = mlm_plot, filename = "~/Downloads/mlm_results.png", width = 10, height = 5, units = "in")




