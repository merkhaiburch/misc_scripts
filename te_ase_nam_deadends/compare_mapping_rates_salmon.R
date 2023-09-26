# Plot distribution of mapping rate between different methods
library(dplyr)
library(ggplot2)

# Set directory
setwd("~/../mbb262-admin/Library/CloudStorage/Box-Box/git_projects/te_ase_nam")

merged_decoy_no500 <- read.delim("data/temp_storage/mapping_rate_salmon_iterations/multiqc_salmon_mergedno500.txt") %>% select(Sample, percent_mapped)
paired_decoy_no500 <- read.delim("data/temp_storage/mapping_rate_salmon_iterations/multiqc_salmon_pairedno500.txt") %>% select(Sample, percent_mapped)
paired_decoy_500buffer <- read.delim("data/temp_storage/mapping_rate_salmon_iterations/multiqc_salmon_paired500.txt") %>% select(Sample, percent_mapped)

# Add metadata
merged_decoy_no500$id <- rep("Merged", nrow(merged_decoy_no500))
paired_decoy_no500$id <- rep("Paired", nrow(paired_decoy_no500))
paired_decoy_500buffer$id <- rep("Paired + 500bp", nrow(paired_decoy_500buffer))

# Combine
together <- rbind(merged_decoy_no500, paired_decoy_no500, paired_decoy_500buffer)

# Get means by group
plot_means <- together %>% group_by(id) %>% summarise(mean = mean(percent_mapped))

# Plot
lala <- ggplot(together, aes(x = percent_mapped, fill = id)) +
  geom_histogram(position="identity", alpha=0.5) +
  geom_vline(data=plot_means, aes(xintercept=mean, color=id), linetype="solid") +
  xlab("Percent Mapped") +
  ylab("Count")+
  theme(legend.position = "bottom")

ggsave(plot = lala, filename = "figs/compare_mapping_rates.png")
