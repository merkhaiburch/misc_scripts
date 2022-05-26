# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-02-23 
# Updated... 2022-05-26

# Description 
# Plot Manhattan-style plot of counts of pleiotropy scores
# And create thresholds of 'pleiotropic' for further investigation
# ---------------------------------------------------------------

# ------------------
# Load in packages
# ------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(dplyr)


# -------------------------------
# Data gathering & formatting
# -------------------------------

# Make shortcut for path
data_dir <- "/Volumes/merrittData1/pleiotropy/interval_data/"
# data_dir <- "~/git_projects/pleiotropy/data"

# Load in unmelted data
all_interval_data <- data.table::fread(paste0(data_dir, "/all_interval_data_unique_filter.csv"))

# Get gene descriptions through gff file
gff <- ape::read.gff("~/Downloads/Zea_mays.B73_RefGen_v4.49.gff3.gz",
                     na.strings = c(".", "?"),
                     GFF3 = TRUE)

# Remove extra columns
gff <- gff %>%
  filter(type == "gene", source == "gramene", seqid %in% c(1,2,3,4,5,6,7,8,9,10)) %>%
  select("seqid", "start", "end", "attributes")

# Parse out gene id from attributes column
gff$gene <- gsub("ID=gene:", "", gff$attributes)
gff$gene <- gsub(";.*", "", gff$gene)

# parse out descriptions
gff$description <- gsub(".*description=", "", gff$attributes)
gff$description <- gsub(";.*", "", gff$description)

# Plotting variables ----------------------------

axis_text_size <- 15
title_text_size <- 16
tag_size <- 16
legend_text_size <- 14
legend_shape_size <- 0.6

# -----------------------------------------------------------------------------------------------------------
# Make Manhattan-style plots - proportion data
# -----------------------------------------------------------------------------------------------------------

# Collect ratios and genes
pleio_ratios <- all_interval_data[,c(1:4,11,38,39,68,70, 81, 92,103:109,145)]

# NAM field --------------------------------------------------------------------------------
temp <- all_interval_data %>% select("rr_id", "start", "seqid", "v4_gene","range_type", 
                                     "nam_filtered_all_true_uniqueCount", 
                                     "npp_unique_perm_1","npp_unique_perm_2","npp_unique_perm_3",
                                     "npp_unique_perm_4","npp_unique_perm_5","npp_unique_perm_6",
                                     "npp_unique_perm_7","npp_unique_perm_8","npp_unique_perm_9",
                                     "npp_unique_perm_10")
dim(temp)

# Take mean over all permutations
temp$perm_mean <- rowMeans(temp[,7:16])

# Calculate the log(observed/permuted mean)
temp$scaled_pleio_mean <- log((temp[,6] + 1)/(temp$perm_mean + 1))
pleio_ratios <- merge(pleio_ratios, temp[,c(1,18)], by = "rr_id")
colnames(pleio_ratios)[20] <- paste0(colnames(pleio_ratios)[20], "_nam_field")

# plot
no <- ggplot(temp, aes(x = start, y = scaled_pleio_mean, color = as.factor(seqid)))+
  geom_point() +
  geom_hline(yintercept =quantile(temp$scaled_pleio_mean, .99)) +
  theme_classic() +
  facet_grid(.~seqid, scales = "free") +
  xlab("Position") +
  ylab("log((Observed + 1)/(Mean Permuted Trait Count + 1))")+
  ggtitle("NAM Field") +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size),
        axis.text.x = element_text(angle = 45, hjust=1),
        plot.tag=element_text(size=tag_size),
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm')) +
  labs(tag = "A") 


# Goodman field  ------------------------------------------------------------------
temp <- all_interval_data %>% select("rr_id", "start", "seqid", "v4_gene", "range_type", "goodman_filtered_physiological_true_uniqueCount", 
                                     "gpp_unique_perm_1","gpp_unique_perm_2","gpp_unique_perm_3",
                                     "gpp_unique_perm_4","gpp_unique_perm_5","gpp_unique_perm_6",
                                     "gpp_unique_perm_7","gpp_unique_perm_8","gpp_unique_perm_9",
                                     "gpp_unique_perm_10")
dim(temp)

# Take mean over all permutations
temp$perm_mean <- rowMeans(temp[,7:16])

# Calculate the log(observed/permuted mean)
temp$scaled_pleio_mean <- log((temp[,6] + 1)/(temp$perm_mean + 1))
pleio_ratios <- merge(pleio_ratios, temp[,c(1,18)], by = "rr_id")
colnames(pleio_ratios)[21] <- paste0(colnames(pleio_ratios)[21], "_goodman_field")

gfo <- ggplot(temp, aes(x = start, y = scaled_pleio_mean, color = as.factor(seqid)))+
  geom_point() +
  geom_hline(yintercept =quantile(temp$scaled_pleio_mean, .99)) +
  theme_classic() +
  facet_grid(.~seqid, scales = "free") +
  xlab("Position") +
  ylab("log((Observed + 1)/(Mean Permuted Trait Count + 1))")+
  ggtitle("Goodman Field") +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size),
        axis.text.x = element_text(angle = 45, hjust=1),
        plot.tag=element_text(size=tag_size),
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm')) +
  labs(tag = "B") 


# Goodman mass features ------------------------------------------------------------------
temp <- all_interval_data %>% select("rr_id", "start", "seqid", "v4_gene", "range_type", 
                                     "goodman_filtered_metabolite_true_uniqueCount", 
                                     "gpm_perm_1","gpm_perm_2","gpm_perm_3",
                                     "gpm_perm_4","gpm_perm_5","gpm_perm_6",
                                     "gpm_perm_7","gpm_perm_8","gpm_perm_9",
                                     "gpm_perm_10")
dim(temp)

# Take mean over all permutations
temp$perm_mean <- rowMeans(temp[,7:16])

# Calculate the log(observed/permuted mean)
temp$scaled_pleio_mean <- log((temp[,6] + 1)/(temp$perm_mean + 1))
pleio_ratios <- merge(pleio_ratios, temp[,c(1,18)], by = "rr_id")
colnames(pleio_ratios)[22] <- paste0(colnames(pleio_ratios)[22], "_goodman_mass")

# plot
gmo <- ggplot(temp, aes(x = start, y = scaled_pleio_mean, color = as.factor(seqid)))+
  geom_point() +
  geom_hline(yintercept =quantile(temp$scaled_pleio_mean, .99)) +
  theme_classic() +
  facet_grid(.~seqid, scales = "free") +
  xlab("Position") +
  ylab("log((Observed + 1)/(Mean Permuted Trait Count + 1))")+
  ggtitle("Goodman Mass Features") +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size),
        axis.text.x = element_text(angle = 45, hjust=1),
        plot.tag=element_text(size=tag_size),
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm')) +
  labs(tag = "C") 


# Goodman expression ------------------------------------------------------------------
# Subset columns
temp <- all_interval_data[,c(1,2,3,11,68,103:144,145:150)]
dim(temp)

# Take mean over all permutations for a single tissue
temp <- data.frame(temp)
temp$GRoot_perm_mean <- rowMeans(temp[grepl(".*unique_GRoot.*", colnames(temp))])
temp$GShoot_perm_mean <- rowMeans(temp[grepl(".*unique_GShoot.*", colnames(temp))])
temp$Kern_perm_mean <- rowMeans(temp[grepl(".*unique_Kern.*", colnames(temp))])
temp$L3Base_perm_mean <- rowMeans(temp[grepl(".*unique_L3Base.*", colnames(temp))])
temp$L3Tip_perm_mean <- rowMeans(temp[grepl(".*unique_L3Tip.*", colnames(temp))])
temp$LMAD_perm_mean <- rowMeans(temp[grepl(".*unique_LMAD.*", colnames(temp))])
temp$LMAN_perm_mean <- rowMeans(temp[grepl(".*unique_LMAN.*", colnames(temp))])

# Calculate the log ratio log(observed/permuted mean)
temp$ratio_GRoot <- log((temp$GRoot_counts + 1)/(temp$GRoot_perm_mean + 1))
temp$ratio_GShoot <- log((temp$GShoot_counts + 1)/(temp$GShoot_perm_mean + 1))
temp$ratio_Kern <- log((temp$Kern_counts + 1)/(temp$Kern_perm_mean + 1))
temp$ratio_L3Base <- log((temp$L3Base_counts + 1)/(temp$L3Base_perm_mean + 1))
temp$ratio_L3Tip <- log((temp$L3Tip_counts + 1)/(temp$L3Tip_perm_mean + 1))
temp$ratio_LMAD <- log((temp$LMAD_counts + 1)/(temp$LMAD_perm_mean + 1))
temp$ratio_LMAN <- log((temp$LMAN_counts + 1)/(temp$LMAN_perm_mean + 1))

# do for count across all tissues
temp$perm_mean <- rowMeans(temp[,49:53])

# Calculate the log(observed/permuted mean)
temp$scaled_pleio_mean <- log((temp$goodman_filtered_expression_uniqueCount + 1)/(temp$perm_mean + 1))

# Just collect the necessary colummns with ratios
temp2 <- temp %>% select("rr_id", "start", "seqid", "scaled_pleio_mean", "ratio_GRoot", "ratio_GShoot", "ratio_Kern",
                         "ratio_L3Base", "ratio_L3Tip", "ratio_LMAD", "ratio_LMAN")
pleio_ratios <- merge(pleio_ratios, temp2[,-c(2:3)], by = "rr_id")
colnames(pleio_ratios)[23:30] <- paste0(colnames(pleio_ratios)[23:30], "_goodman_expression")

# melt data
temp2 <- tidyr::pivot_longer(temp2,
                             cols = ratio_GRoot:ratio_LMAN,
                             names_to = "tissue_type",
                             values_to = "ratio")

# plot
geo <- ggplot(temp2, aes(x = start, y = ratio, color = as.factor(seqid)))+
  geom_point() +
  geom_hline(yintercept =quantile(temp2$ratio, .99)) +
  theme_classic() +
  facet_grid(.~seqid, scales = "free") +
  xlab("Position") +
  ylab("log((Observed + 1)/(Mean Permuted Trait Count + 1))")+
  ggtitle("Goodman Expression") +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size),
        axis.text.x = element_text(angle = 45, hjust=1),
        plot.tag=element_text(size=tag_size),
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm')) +
  labs(tag = "D") 


# Combine all plots & export------------------------------------------------------------------------
ggsave("~/git_projects/pleiotropy/images/manhattan_plots_pleiotropy_scores.png",
       plot = ggpubr::ggarrange(no, gfo, gmo, geo,
                                nrow = 2, ncol = 2,
                                common.legend = TRUE, legend = "bottom", align = "v"),
       width = 20,
       height = 15,
       units = "in",
       dpi = "retina")

