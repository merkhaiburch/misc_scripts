# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-10-27 
#
# Description 
#   - Checking flowering time triplet haplotye peak calls for 
#   - rentention and accuracy in the NAM and Goodman Panels
# ---------------------------------------------------------------

# read in required packages
library(data.table)
library(patchwork)
library(dplyr)

#create a list of the files from your target directory
nam_files <- list.files("~/Documents/mgwas_results/ft_hap_peaks/nam/")
goodman_files <- list.files("~/Documents/mgwas_results/ft_hap_peaks/goodman/")
nam_files_all <- list.files("~/Documents/mgwas_results/ft_hap_peaks/unfiltered_nam/")
goodman_files_all <- list.files("~/Documents/mgwas_results/ft_hap_peaks/unfiltered_goodman/")


# initiate a blank data frame
nam_dat <- data.frame()
goodman_dat <- data.frame()
nam_dat_all <- data.frame()
goodman_dat_all <- data.frame()

# Rbind all chromosomal results
setwd("~/Documents/mgwas_results/ft_hap_peaks/nam")
for (i in 1:length(nam_files)){
  temp_data <- data.table::fread(nam_files[i], stringsAsFactors = F) 
  nam_dat <- rbindlist(list(nam_dat, temp_data), use.names = T)
}

setwd("~/Documents/mgwas_results/ft_hap_peaks/goodman/")
for (i in 1:length(goodman_files)){
  temp_data <- data.table::fread(goodman_files[i], stringsAsFactors = F) 
  goodman_dat <- rbindlist(list(goodman_dat, temp_data), use.names = T)
}

setwd("~/Documents/mgwas_results/ft_hap_peaks/unfiltered_nam")
for (i in 1:length(nam_files_all)){
  temp_data <- data.table::fread(nam_files_all[i], stringsAsFactors = F) 
  nam_dat_all <- rbindlist(list(nam_dat_all, temp_data), use.names = T)
}

setwd("~/Documents/mgwas_results/ft_hap_peaks/unfiltered_goodman")
for (i in 1:length(goodman_files_all)){
  temp_data <- data.table::fread(goodman_files_all[i], stringsAsFactors = F) 
  goodman_dat_all <- rbindlist(list(goodman_dat_all, temp_data), use.names = T)
}

# Change column names of results containing all p-values
colnames(nam_dat_all)[8:9] <- c("snp_coord", "seqid") 
colnames(goodman_dat_all)[3:4] <- c("seqid","snp_coord") 

# Plot results NAM
nam_dat <- nam_dat %>% filter(Trait == "DTS_BLUP.Buckler_2009")
nam_dat_all <- nam_dat_all %>% filter(Trait == "DTS_BLUP.Buckler_2009")
a1  <- manhattan_without_windows8(results = nam_dat, upper_bound = 30,
                                  title = "NAM: Filtered top 1% SNPs in sliding triplet haplotype \n Chrom 8: DTS ~ SNP + 3 PCs")
a2 <- manhattan_without_windows10(results = nam_dat, upper_bound = 70,
                                  title = "NAM: Filtered top 1% SNPs in sliding triplet haplotype \n Chrom 10: DTS ~ SNP + 3 PCs") 

b1  <- manhattan_without_windows8(results = nam_dat_all, upper_bound = 30,
                                  title = "NAM: All SNP results \n Chrom 8: DTS ~ SNP + 3 PCs")
b2 <- manhattan_without_windows10(results = nam_dat_all, upper_bound = 70,
                                  title = "NAM: All SNP results \n Chrom 10: DTS ~ SNP + 3 PCs ") 

# Plot all at once
# (a1/a2)|(b1/b2)
ggsave("~/Documents/mgwas_results/ft_hap_peaks/nam_filtered_unfiltered_results.png", 
       (a1/a2)|(b1/b2))

# Plot results Goodman
goodman_dat <- goodman_dat %>% filter(Trait == "DTS_BLUP.Buckler_2009")
goodman_dat_all <- goodman_dat_all %>% filter(Trait == "DTS_BLUP.Buckler_2009")

g1  <- manhattan_without_windows8(results = goodman_dat, upper_bound = 25,
                                  title = "Goodman: Filtered top 1% SNPs in sliding triplet haplotype \n Chrom 8: DTS ~ SNP + 3 PCs")
g2 <- manhattan_without_windows10(results = goodman_dat, upper_bound = 25,
                                  title = "Goodman: Filtered top 1% SNPs in sliding triplet haplotype \n Chrom 10: DTS ~ SNP + 3 PCs") 

g3  <- manhattan_without_windows8(results = goodman_dat_all, upper_bound = 25,
                                  title = "Goodman: All SNP results \n Chrom 8: DTS ~ SNP + 3 PCs")
g4 <- manhattan_without_windows10(results = goodman_dat_all, upper_bound = 25,
                                  title = "Goodman: All SNP results \n Chrom 10: DTS ~ SNP + 3 PCs ") 

# Plot all at once
# (g1/g2)|(g3/g4)
ggsave("~/Documents/mgwas_results/ft_hap_peaks/goodman_filtered_unfiltered_results.png", 
       (g1/g2)|(g3/g4))


# Functions to plot
manhattan_without_windows10 <- function(results, title, upper_bound){
  
  # just chromosome 10
  chr_results <- results %>% filter(seqid == 10)
  
  # Get significance thresholds
  sig_threshold_1 = -log10(quantile(results$p, probs = 0.05))
  sig_threshold_2 = 8
  
  # Text position
  text_pos <- upper_bound - 4
  
  # Plot the thing
  ggplot(chr_results, aes(x = snp_coord, y = -log10(p), size = -log10(p))) +
    geom_point(size = 1) +
    geom_hline(yintercept = sig_threshold_1, color = "grey40", linetype = "dashed") +
    geom_hline(yintercept = sig_threshold_2, color = "black", linetype = "dashed") +
    scale_y_continuous(expand = c(0,0), limits = c(1.5, upper_bound)) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, y = "-log10(p)", title = title) + 
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 9, vjust = 0.5),
      axis.text.y = element_text(size = 9, vjust = 0.5)) +
    annotate("text", x = 77700138, y = text_pos, label = "ZmCCA1") +
    geom_vline(xintercept=77700138) +
    annotate("text", x = 77722777, y = text_pos-2, label = "ZmLHY1") +
    geom_vline(xintercept=77722777) +
    annotate("text", x =94430850 , y = text_pos, label = "ZmCCT") +
    geom_vline(xintercept=94430850) +
    annotate("text", x = 141561862, y = text_pos, label = "ZFL1") +
    geom_vline(xintercept=141561862)
  
}
manhattan_without_windows8 <- function(results, title, upper_bound){
  
  # Filter results to just chrom 8
  chr_results <- results %>% filter(seqid == 8)
  
  # Get significance thresholds
  sig_threshold_1 = -log10(quantile(results$p, probs = 0.05))
  sig_threshold_2 = -log10(quantile(results$p, probs = 0.01))
  
  # Text position
  text_pos <- upper_bound - 4
  
  # Plot the thing
  ggplot(chr_results, aes(x = snp_coord, y = -log10(p), size = -log10(p))) +
    geom_point(size = 1) +
    geom_hline(yintercept = sig_threshold_1, color = "grey40", linetype = "dashed") +
    geom_hline(yintercept = sig_threshold_2, color = "black", linetype = "dashed") +
    scale_y_continuous(expand = c(0,0), limits = c(1.5, upper_bound)) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, y = "-log10(p)", title = title) + 
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 9, vjust = 0.5),
      axis.text.y = element_text(size = 9, vjust = 0.5)) +
    annotate("text", x = 21787238, y = text_pos, label = "GIGZ1A") +
    geom_vline(xintercept=21787238) +
    annotate("text", x = 126880531, y = text_pos, label = "ZCN8") +
    geom_vline(xintercept=126880531) +
    annotate("text", x = 135946643, y = text_pos, label = "vgt1") +
    geom_vline(xintercept=135946643) +
    annotate("text", x =136009216 , y = text_pos-2, label = "ZmRAP2.7") +
    geom_vline(xintercept=136009216) +
    annotate("text", x = 163614240, y = text_pos, label = "ZmHy2") +
    geom_vline(xintercept=163614240)
}
