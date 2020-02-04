# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-02-03 
#
# Description 
#   - Test effect of collapsing PHG haplotypes with 26 NAM founders
#   - maxDiv levels = 0.0001, 0.0005, 0.005, 0.001, 0.01
# ---------------------------------------------------------------

# Set workdir
setwd("~/Box Sync/Cornell_PhD/labProjects/hackathons/2020_02_03_hackathon")

# Libraries
library(ggplot2)
library(dplyr)

# Load file
# Had to delete 3,5, and some empty columns to the right for file to read in
# In sublime I removed by ^\s*$, ^\s, and by DarrowSomething
# hap <- read.table("haplotype_table_clean.txt", header = F, sep=" ") # Evan generated

# Add header
# colnames(hap) <- c("haplotype_id", "gamete_grp_id", "ref_range_id", "seq_len")

# Plotting by seq_len
hist(hap$seq_len)


# ------------------
# Evan's dataframe
# ------------------

# Taxa ID by haplotype seqeunce length
ggplot(hap, aes(x = as.factor(gamete_grp_id), y = seq_len)) + 
  geom_boxplot()

# Histogram of all seq_lens
hist(hap$seq_len)

# Genic haplotypes only
genic <- hap[hap$ref_range_id < 35677, ]
ggplot(genic, aes(x = as.factor(gamete_grp_id), y = seq_len)) + 
  geom_boxplot() +
  ylim(0,50000)

# Number of haplotypes per reference range and other basic stats by ref range
temp <- genic %>%
  group_by(ref_range_id) %>%
  summarise(mean = mean(seq_len), min = min(seq_len), max = max(seq_len), n = n())


# -----------------------
# Uncollapsed haplotypes
# -----------------------

# No modification needed
hap <- read.table("cintaMap.txt", header = T) # Lynn generated

# Melt data 
library(reshape)
hap$annotation <- c(rep("genic", 71354/2), rep("intergenic", 71354/2)) # annotation
mdata <- melt(hap, id=c("refRangeID","chr","start","end", "annotation"))

# Count number of missing haplotypes per taxa
temp <- colSums(hap == 0)

# Plot faceted lengths  of haplotypes by taxa
ggplot(mdata, aes(x = variable, y = log10(value), fill = annotation)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_blank(),
        text = element_text(size=20)) +
  labs(x="Taxa", y = "log10(haplotype length)") +
  facet_wrap(~annotation, nrow = 2)+ guides(fill=FALSE)

