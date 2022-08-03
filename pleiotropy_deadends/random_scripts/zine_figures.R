# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-05-24
# Updated... 2022-05-30
# 
# Description 
# Plotting unadjusted count of the number of uniue traits per interval
# Plotting:
# ]Filtered vs permuted unadjusted data with genic and intergenic
# ---------------------------------------------------------------

# load packages
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggbreak)

# Make shortcut for path
data_dir <- "/Volumes/merrittData1/pleiotropy/interval_data"

# Plotting variables ----------------------------

axis_text_size <- 4.5
title_text_size <- 5
tag_size <- 4.5

# The color blind palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# ---------------------------------------------------------------------------
# 1) Plot un-adjusted, splitting by genic or intergenic intervals
# log(observed data + 1/mean permuted data count + 1)
# ---------------------------------------------------------------------------
# all_interval_data <- data.table::fread(paste0(data_dir, "/all_interval_data_unique_filter.csv"))

# Physiological NAM-------------------------------------------------
temp <- all_interval_data %>% select("rr_id", "range_type", "nam_filtered_all_true_uniqueCount", 
                                     "npp_unique_perm_1","npp_unique_perm_2","npp_unique_perm_3",
                                     "npp_unique_perm_4","npp_unique_perm_5","npp_unique_perm_6",
                                     "npp_unique_perm_7","npp_unique_perm_8","npp_unique_perm_9",
                                     "npp_unique_perm_10")
dim(temp)

# Take mean over all permutations
temp$perm_mean <- rowMeans(temp[,4:13])

# Calculate the log(observed/permuted mean)
temp$ratio_pleio_mean <- (temp[,3] + 1)/(temp$perm_mean + 1)

# Calculate the log(observed/permuted mean)
temp$scaled_pleio_mean <- log((temp$nam_filtered_all_true_uniqueCount + 1)/(temp$perm_mean + 1))

# check stats
temp %>% group_by(range_type) %>% summarise(m = mean(scaled_pleio_mean), med = median(scaled_pleio_mean), n = n())

# plot
pn <- ggplot(data=temp, aes(x = scaled_pleio_mean, fill = range_type)) +
  geom_histogram(size=0.5, bins = 100, alpha = 0.4) +
  xlab("log(Observed/Permuted)") +
  ylab("Count of Intervals")+
  ggtitle("NAM Field") +
  theme_classic()+
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="none", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size), 
        axis.ticks = element_line(colour = "black", size = 0.3),
        axis.line = element_line(colour = "black", size = 0.3)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  labs(tag = "A") +
  ylim(0, 17000) +
  geom_vline(xintercept = log(6/2), linetype = "dashed", size = 0.3) 
 


# Physiological Goodman---------------------------------------------
temp <- all_interval_data %>% select("rr_id", "range_type", "goodman_filtered_physiological_true_uniqueCount", 
                                     "gpp_unique_perm_1","gpp_unique_perm_2","gpp_unique_perm_3",
                                     "gpp_unique_perm_4","gpp_unique_perm_5","gpp_unique_perm_6",
                                     "gpp_unique_perm_7","gpp_unique_perm_8","gpp_unique_perm_9",
                                     "gpp_unique_perm_10")
dim(temp)

# Take mean over all permutations
temp$perm_mean <- rowMeans(temp[,4:13])

# Calculate the log(observed/permuted mean)
temp$scaled_pleio_mean <- log((temp[,3] + 1)/(temp$perm_mean + 1))

# check stats
temp %>% group_by(range_type) %>% summarise(m = mean(scaled_pleio_mean), med = median(scaled_pleio_mean), n = n())

# plot
pg <- ggplot(data=temp, aes(x = scaled_pleio_mean, fill = range_type)) +
  geom_histogram(size=0.5, bins = 100, alpha = 0.4) +
  xlab("log(Observed/Permuted)") +
  ylab("Count of Intervals")+
  ggtitle("Goodman Field") +
  theme_classic()+
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="none", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size), 
        axis.ticks = element_line(colour = "black", size = 0.3),
        axis.line = element_line(colour = "black", size = 0.3)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  labs(tag = "B") +
  ylim(0, 17000) +
  geom_vline(xintercept = log(6/2), linetype = "dashed", size = 0.3) 

# Metabolite-------------------------------------------------------
temp <- all_interval_data %>% select("rr_id", "range_type", "goodman_filtered_metabolite_true_uniqueCount", 
                                     "gpm_perm_1","gpm_perm_2","gpm_perm_3",
                                     "gpm_perm_4","gpm_perm_5","gpm_perm_6",
                                     "gpm_perm_7","gpm_perm_8","gpm_perm_9",
                                     "gpm_perm_10")
dim(temp)

# Take mean over all permutations
temp$perm_mean <- rowMeans(temp[,4:13])

# Calculate the log(observed/permuted mean)
temp$scaled_pleio_mean <- log((temp[,3] + 1)/(temp$perm_mean + 1))

# check stats
temp %>% group_by(range_type) %>% summarise(m = mean(scaled_pleio_mean), med = median(scaled_pleio_mean), n = n())

# plot
m <- ggplot(data=temp, aes(x = scaled_pleio_mean, fill = range_type)) +
  geom_histogram(size=0.5, bins = 100, alpha = 0.4) +
  xlab("log(Observed/Permuted)") +
  ylab("Count of Intervals")+
  ggtitle("Goodman Mass Features") +
  theme_classic()+
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="none", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size), 
        axis.ticks = element_line(colour = "black", size = 0.3),
        axis.line = element_line(colour = "black", size = 0.3)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  labs(tag = "C") +
  ylim(0, 17000) +
  geom_vline(xintercept = log(6/2), linetype = "dashed", size = 0.3) 

# Expression--------------------------------------------------------
# Subset columns
temp <- all_interval_data[,c(1,68,103:144)]
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

# Just collect the necessary colummns with ratios
temp2 <- temp %>% select("rr_id", "range_type", "ratio_GRoot", "ratio_GShoot", "ratio_Kern",
                         "ratio_L3Base", "ratio_L3Tip", "ratio_LMAD", "ratio_LMAN")

# melt data
temp2 <- tidyr::pivot_longer(temp2,
                             cols = ratio_GRoot:ratio_LMAN,
                             names_to = "tissue_type",
                             values_to = "ratio")

# plot
e <- ggplot(data=temp2, aes(x = ratio, fill = range_type)) +
  geom_histogram(size=0.5, bins = 100, alpha = 0.4) +
  xlab("log(Observed/Permuted)") +
  ylab("Count of Intervals")+
  ggtitle("Goodman Expression") +
  theme_classic()+
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="none", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size), 
        axis.ticks = element_line(colour = "black", size = 0.3),
        axis.line = element_line(colour = "black", size = 0.3)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  labs(tag = "D") +
  ylim(0, 85000) +
  geom_vline(xintercept = log(6/2), linetype = "dashed", size = 0.3) 

# Save to file-------------------------------------------------------
ggsave(filename = "~/git_projects/pleiotropy/images/zine_proportion_hits_original_vs_mean_permuted_log.png",
       plot = ggpubr::ggarrange(pn, pg, m, e,
                                nrow = 2, ncol = 2, 
                                common.legend = TRUE, 
                                legend = "none", 
                                align = "v"),
       width = 2.7,
       height = 2.2,
       units = "in",
       dpi = "retina"
)


