# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-11-05
#
# Description 
#   - Look at extent of pleiotropy from all SNPs and and triplet filtered SNPs
#   - compare results with physiological traits 
# Sister script: meta_results_over_regions.R
# ---------------------------------------------------------------

# Load packages
library(data.table)
library(ggplot2)
library(logging)
library(magrittr)
library(patchwork)
library(tibble)
library(tidyr)
library(dplyr)


## 'Hack' for getting scientific notation
addUnits <- function(n) {
  labels <- ifelse(n < 1000, n,  # less than thousands
                   ifelse(n < 1e6, paste0(round(n/1e3), 'k'),  # in thousands
                          ifelse(n < 1e9, paste0(round(n/1e6), 'M'),  # in millions
                                 ifelse(n < 1e12, paste0(round(n/1e9), 'B'), # in billions
                                        ifelse(n < 1e15, paste0(round(n/1e12), 'T'), # in trillions
                                               'too big!'
                                        )))))
  return(labels)
}


# ---------------------------------------------
# Get data from the following sources from CBSU
# ---------------------------------------------

# Data in different genic regions and genome wide for Goodman panel eQTL is in
# /data1/users/bm646/askdb_output_data_kremling_eqtl

# NAM and Goodman Panel data is here:


# -----------------------
# Load & Process the data
# -----------------------

# eQTL Triplet filterd results with >0 traits mapping to single SNP
trip_files <- list.files("/workdir/mbb262/result_counts/", pattern = "count_data_all_*")
trip_short <- gsub("count_data_all_", "", trip_files)
trip_short <- gsub(".csv", "", trip_short)
setwd("/workdir/mbb262/result_counts/")
filt_dat <- data.frame() # Empty df
for (i in 1:length(trip_files)){
  temp_data <- data.table::fread(trip_files[i], stringsAsFactors = F) 
  temp_data$tissue <- rep(trip_short[i], nrow(temp_data))
  filt_dat <- rbindlist(list(filt_dat, temp_data), use.names = T)
}

## Load NAM and Goodman data names
pleio_files_pop <- list.files(path = "local_data/", pattern = ".csv", full.names = TRUE)
pleio_files_nam <- list.files(path = "local_data/hapmap_coords_nam/", full.names = TRUE)

## Go through each chromosome (not robust) ----
goodman_ls <- list()
for (i in 1:10) {
  logging::loginfo(paste0("Parsing: ", pleio_files_pop[i]))
  
  # Load
  pleio_pop <- data.table::fread(pleio_files_pop[i])
  pleio_nam <- data.table::fread(pleio_files_nam[i])
  
  # JOIN
  data.table::setkey(pleio_pop, snp_coord)
  data.table::setkey(pleio_nam, snp_coord)
  goodman_vars <- pleio_pop[pleio_nam, nomatch = 0]
  
  # Clean up columns
  goodman_vars <- goodman_vars[, -6]
  
  goodman_ls[[i]] <- goodman_vars
}


## Row bind all data from list elements ----
goodman_df <- do.call("rbind", goodman_ls)
data.table::fwrite(goodman_df, file = "local_data/goodman_df.csv")


## Convert to tibble (for Tidyverse stuff) ----
goodman_tib <- tibble::as_tibble(goodman_df)
goodman_tib <- goodman_tib[, -3]


## Pivot data ----
goodman_pivot <- goodman_tib %>%
  pivot_longer(!c(snp_seqid, snp_coord), "pop_id")


# ----------------------------------------------
# Are the same set of markers consistently the 
# outliers/have the most counts?
# ----------------------------------------------

# filt_dat has all tissues

# Get top 1000 hits (i.e. markers with the most trait associations) for each tissue
top_1k <- filt_dat %>% group_by(tissue) %>% slice_max(order_by = count, n = 1000)

# Try sorting with datatable for speed
# library(data.table)
# temp <- filt_dat[order(-rank(count))]


# Make each tissue it's own element within a list, this is hardcoded and not the best
marker_list <- list(top_1k %>% filter(tissue == "GRoot") %>% pull(Marker),
                    top_1k %>% filter(tissue == "GShoot") %>% pull(Marker),
                    top_1k %>% filter(tissue == "Kern") %>% pull(Marker),
                    top_1k %>% filter(tissue == "L3Base") %>% pull(Marker),
                    top_1k %>% filter(tissue == "L3Tip") %>% pull(Marker),
                    top_1k %>% filter(tissue == "LMAD") %>% pull(Marker),
                    top_1k %>% filter(tissue == "LMAN") %>% pull(Marker))
names(marker_list) <- c("GRoot", "GShoot","Kern","L3Base", "L3Tip","LMAD","LMAN")

temp <- marker_list[-1]
# Plot overlap between tissues
library(VennDiagram)
venn.plot <- venn.diagram(
  x = marker_list,
  filename = "/workdir/mbb262/test.png"
)

# Looks like many markers are repeated, try upping MAF filter for single tissue, remap, and reassess results

# Try cheap method and remove the top 1000 SNPs with the most traits mapping to them
#   and see what the distributions look like



# -------------------
# Calculate medians
# -------------------

# Calculate median pleiotropy across tissues
tissue_triplet_filt_more0 <- filt_dat %>% 
  group_by(tissue) %>% 
  summarise(count = median(count))

# Extract row with largest value
filt_dat[which.max(filt_dat$count),]

# Sort df by count
temp <- as.data.table(filt_dat) %>% 
  arrange(desc(count)) %>% 
  filter(tissue == "L3Base")

box <- ggplot(temp, aes(x = seqname, y = count)) +
  geom_boxplot() +
  xlab("# eQTLs associated with a single SNP by chromosome for L3Base") +
  ylab("Number of SNPs")

ggsave(
  filename = "pleio_eQTL_L3Base.png",
  plot = box,
  width = 6.5,
  height = 4.5,
  units = "in",
  dpi = 300
)

# -----------------------
# Plot results
# -----------------------

# Plot results
vis_A <- filt_dat %>%
  ggplot() +
  aes(x = count, fill = tissue) +
  geom_histogram(position = "identity", bins = 200, alpha = 0.6) +
  xlab("# eQTLs associated with a single SNP") +
  ylab("Number of SNPs") +
  scale_y_continuous(labels = addUnits) +
  theme(legend.position = "bottom")

ggsave(
  filename = "pleio_all_eQTL.png",
  plot = vis_A,
  width = 6.5,
  height = 4.5,
  units = "in",
  dpi = 300
)


# Plot this with physiological traits?



# -----------------------------------------------------
# Pseudocode for counting the number of traits per SNP
# outside of askdb
# -----------------------------------------------------


