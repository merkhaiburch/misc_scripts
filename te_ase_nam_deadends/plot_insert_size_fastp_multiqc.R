# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-08-09
# Updated... 2023-08-09
#
# Description:
# Plot insert sizes for each sample while dropping outliers
# data from a multiqc output on fastp json files
# ------------------------------------------------------------------------------


library(ggplot2)
library(dplyr)
library(tidyr)

insertSize <- read.delim("data/temp_storage/mqc_fastp-insert-size-plot_1.txt", header = F)

insertSize <- insertSize[grepl("MS21",insertSize[,1]),]

colnames(insertSize) <- c("sample", rep(1:(ncol(insertSize)-1)))
insertSize$sample <- gsub("_.*", "", insertSize$sample)

meltInsertSize <- insertSize %>% 
  tidyr::pivot_longer(cols = !sample, names_to = "insert_size", values_to = "read_percent")
dim(meltInsertSize)

before <- ggplot(meltInsertSize, aes(x = insert_size, y = read_percent)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

lala <- read.csv("data/bad_taxa.txt")
meltInsertSize <- meltInsertSize[!(meltInsertSize$sample %in% lala$short_name),]
meltInsertSize$insert_size <- as.numeric(as.character(meltInsertSize$insert_size))

b <- ggplot(meltInsertSize, aes(x = insert_size, y = read_percent, group = sample)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
b
ggsave("figs/insert_size.png", b)

meltInsertSize %>% filter(insert_size < 50 & read_percent > 2)
