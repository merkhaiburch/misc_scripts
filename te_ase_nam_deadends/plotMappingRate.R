

# Plot histogram of salmon mapping rate for inbreds

# Set working directory
setwd("/Users/mbb262-admin/Library/CloudStorage/Box-Box/git_projects/te_ase_nam/data/temp_storage")

# Load packages
library(ggplot2)
library(dplyr)

# Load in data
inbred1 <- read.delim("mapping_rate_salmon_iterations/multiqc_salmon_01.txt")
inbred2 <- read.delim("mapping_rate_salmon_iterations/multiqc_salmon_02.txt")
inbred3 <- read.delim("mapping_rate_salmon_iterations/multiqc_salmon_03.txt")
inbred4 <- read.delim("mapping_rate_salmon_iterations/multiqc_salmon_04.txt")
inbred5 <- read.delim("mapping_rate_salmon_iterations/multiqc_salmon_05.txt")
inbred6 <- read.delim("mapping_rate_salmon_iterations/multiqc_salmon_06.txt")
inbred7 <- read.delim("mapping_rate_salmon_iterations/multiqc_salmon_07_new.txt")
inbred8 <- read.delim("mapping_rate_salmon_iterations/multiqc_salmon_08.txt")

# Add identifiers
inbred1$id <- rep("1 - Merged mRNA", nrow(inbred1))
inbred2$id <- rep("5 - Paired liftoff all B73 transcripts", nrow(inbred2))
inbred3$id <- rep("2 - Paired mRNA +500bp", nrow(inbred3))
inbred4$id <- rep("3 - Paired mRNA +500bp no decoy", nrow(inbred4))
inbred5$id <- rep("4 - Paired B73 cDNAs", nrow(inbred5))
inbred6$id <- rep("6 - Paired liftoff max B73 reelGene transcripts", nrow(inbred6))
inbred7$id <- rep("7 - Paired liftoff max transcript each genome", nrow(inbred7))
inbred8$id <- rep("8 - Paired liftoff ALL transcripts each genome", nrow(inbred8))

# Test correlations between two sets
# Its =1, I perfectly re-created cDNAs with liftoff (good)
test <- merge(inbred2, inbred5, by = "Sample")
cor(test$percent_mapped.x, test$percent_mapped.y)

# Combine data
together <- rbind(inbred1, inbred2, inbred3, inbred4, inbred5, inbred6, inbred7, inbred8)

# Make consistent names
together$Sample <- gsub("_.*", "", together$Sample)

# Remove bad one just for now
# together <- together %>% filter(id != "Paired liftoff max transcript each genome")

# Get means by group
mu <- together %>% group_by(id) %>% summarise(grp.mean = mean(percent_mapped))
mu

# Plot mapping rate
library(ggridges)
ridge <- ggplot(together, aes(x = percent_mapped, y = id, group = id, fill = id)) +
  stat_binline(bins = 35, scale = 2.2, draw_baseline = FALSE) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(legend.position = 'none')
ridge
ggsave("/Users/mbb262-admin/Library/CloudStorage/Box-Box/git_projects/te_ase_nam/figs/salmon_inbred_to_founder_mapping_rate_8methods_ridgeline.png", ridge)


# Plot decoy rate
together$decoyrare <- together$num_decoy_fragments/together$num_processed
ahhh2 <- ggplot(together, aes(x = decoyrare))+
  geom_histogram(bins =40) +
  geom_vline(xintercept = mean(together$decoyrare)) +
  annotate("text", x = mean(together$decoyrare)+0.003, y = 9, label = paste0("mean = ", round(mean(together$decoyrare), 2))) +
  xlab("Decoy rate (inbreds to founder)") +
  ylab("Count of Reads")
ahhh2
ggsave("/Users/mbb262-admin/Library/CloudStorage/Box-Box/git_projects/te_ase_nam/figs/salmon_inbred_to_founder_decoy_rate.png", ahhh2)

