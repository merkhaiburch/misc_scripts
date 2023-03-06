# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-10-08
# Updated... 2022-10-08
#
# Description:
# Plot model results 
# ------------------------------------------------------------------------------

# Load packages
library(dplyr)
library(ggplot2)
library(data.table)


## Assemble data ---------------------------------------------------------------

# Load in model effect results
effect_files  <- list.files(pattern = '*_effects_*', path = "/workdir/mbb262/nam_hybrid_rnaseq/output/model_results/", full.names = TRUE)
effect_tables <- lapply(effect_files, data.table::fread, header = TRUE)
combined_model_effects <- do.call(rbind , effect_tables)

# Load in model stats results
stats_files  <- list.files(pattern = '*_stats_*', path = "/workdir/mbb262/nam_hybrid_rnaseq/output/model_results/", full.names = TRUE)
stats_tables <- lapply(stats_files, data.table::fread, header = TRUE)
combined_model_stats <- do.call(rbind , stats_tables)
combined_model_stats <- combined_model_stats %>% filter(Marker != "None")

# Check dimensions before
dim(combined_model_effects)
dim(combined_model_stats)

# Get unique rows  b/c I split this across machines
combined_model_effects <- combined_model_effects %>% unique()
combined_model_stats <- combined_model_stats %>% unique()

# Check dimensions after
dim(combined_model_effects)
dim(combined_model_stats)

# Turn into numeric
combined_model_stats$p <- as.numeric(as.character(combined_model_stats$p))


## Plot results ----------------------------------------------------------------

# Color palette and sizing
crPal <- croix::croix_palette(name = "mov_edward_scissorhands", n = 5)
axis_text_size <- 19

# Make plot
plot1 <- ggplot(combined_model_stats, aes(x = p, fill = tissue)) +
  geom_histogram(position="identity", alpha=0.5, bins = 100) +
  scale_colour_manual(values=crPal) +
  scale_fill_manual(values=crPal) +
  labs(x = "Mixed Model P-value",
       y = "Count") +
  theme_bw() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size))  +
  theme(legend.position="none",
        legend.key.size = unit(1.5, 'cm'),
        legend.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size))+
  facet_wrap(~tissue)
plot1


# Save
ggsave("/home/mbb262/git_projects/te_ase_nam/figs/nam_inbred_expression_model.png", 
       plot1, height = 10, width = 11.5, units = "in")


## Plotting/Counting -----------------------------------------------------------

# Count the number of significant genes across models
combined_model_stats %>%
  filter(p <= 0.05) %>% 
  group_by(tissue) %>% 
  summarize(unique_genes = n_distinct(Trait), unique_tes = n_distinct(Marker))

# Count the number of significant genes across models - across tissues
combined_model_stats %>%
  filter(p <= 0.05) %>% 
  summarize(unique_genes = n_distinct(Trait), unique_tes = n_distinct(Marker))


## Plot effects ----------------------------------------------------------------

# Combine effects and stats
stats_effects <- merge(combined_model_stats, combined_model_effects, 
                       by = c("Trait", "Marker", "tissue"))

# Subset to just significant markers
sub_stats_effects <- stats_effects %>% filter(p <= 0.05)

plot2 <- ggplot(sub_stats_effects, aes(x = p, y = Effect, color = tissue)) +
  geom_point() +
  geom_smooth(method=lm) +
  scale_colour_manual(values=crPal) +
  scale_fill_manual(values=crPal) +
  labs(x = "Mixed Model P-value",
       y = "Effect Size") +
  theme_bw() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size))  +
  theme(legend.position="none",
        legend.key.size = unit(1.5, 'cm'),
        legend.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size))+
  facet_wrap(~tissue)

plot2


# Save
ggsave("/home/mbb262/git_projects/te_ase_nam/figs/nam_inbred_expression_model_effects.png", 
       plot2, height = 10, width = 11.5, units = "in")


