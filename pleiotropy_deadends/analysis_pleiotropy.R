# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-03-30
# Date Updated: 2021-08-03
#
# Description 
#   - Using pleiotropy scores from Shannon diversity index
#   - and doing different analyses
# ---------------------------------------------------------------

# ----------------------------
# Load in packages
# -----------------------------

library(dplyr)
library(data.table)
library(dtplyr)
# library(entropy)
library(ggplot2)
library(lmerTest)
library(lme4)
library(patchwork)
library(tidyr)


# -----------------------------------------------------------------
#                         Gather data
# -----------------------------------------------------------------

# # ----------------------------------
# # Gather genic/intergenic intervals
# # ----------------------------------
# 
# # Get intervals from blfs1
# # scp mbb262@cbsublfs1.tc.cornell.edu:/data1/users/mbb262/haplotype_ranges/genic_intergenic_intervals_b73_v4.49.csv /workdir/mbb262
# 
# # Genic and intergenic ranges (not from the phg but from genic and intergenic ranges in v4 b73)
# ranges <- data.table::fread("/workdir/mbb262/genic_intergenic_intervals_b73_v4.49.csv", header = TRUE)
# 
# # Change column names
# colnames(ranges) <- c("seqid", "start", "end", "rr_id")


# --------------
# Read in gwa counts and entropy data
# --------------

# Location on blfs1
# /data1/users/mbb262/results/pleiotropy_results/diversity_pleiotropy_2021_07/aggregated_by_pop_or_trait/goodman_nam_allTraitTypes_entropy_gwaCounts.txt

# local data directory
data_dir <- "/workdir/mbb262/interval_data/"
# data_dir <- "~/Box/Cornell_PhD/labProjects/hap_gwa/data/interval_data/"
# data_dir <- "C:/Users/merri/Box/Cornell_PhD/labProjects/hap_gwa/data/interval_data/"

# Dataset includes physiological, metabolite, and expression entropy 
#     and the number of GWAS hits per interval
entropy_gwasCounts <- read.csv(paste0(data_dir, "goodman_nam_allTraitTypes_entropy_gwaCounts.txt"))

# Change column names
colnames(entropy_gwasCounts)[1] <- "rr_id"

# Get expression data
expression_int <- data.table::fread(paste0(data_dir, "walley_rna_protein_ranges_overlaps.csv"))

# get LD data
ld_nam <- data.table::fread(paste0(data_dir, "aggregated_ld_nam.csv"), drop = c("seqid", "start", "end"))
colnames(ld_nam)[2] <- "average_r2_ld_nam" 
ld_282 <- data.table::fread(paste0(data_dir, "aggregated_ld_282.csv"), drop = c("seqid", "start", "end"))
colnames(ld_282)[2] <- "average_r2_ld_282" 

# get gerp data
gerp <- data.table::fread(paste0(data_dir, "gerp_aggregated_all_thresholds.csv"))


# Combine different datasets
all_interval_data <- Reduce(function(x, y) merge(x, y, by = "rr_id"), 
                            list(entropy_gwasCounts, expression_int, ld_nam, ld_282, gerp))

# do data cheks
colnames(all_interval_data)
dim(all_interval_data)
head(all_interval_data)

# Add identifier for genic or intergenic ranges
all_interval_data$range_type <- gsub("_rr.*", "", all_interval_data$rr_id)

# calculate interval size
all_interval_data$interval_size <- all_interval_data$end - all_interval_data$start

# rearrange data
all_interval_data <- all_interval_data %>%  
  select(rr_id, range_type, seqid, start, end, interval_size, 
         v4_gene_protein, v4_start_protein, v4_end_protein,
         v4_gene_rna, v4_start_rna, v4_end_rna,
         average_r2_ld_nam, average_r2_ld_282, everything())

# log transform the response, fixed, and random effects
all_interval_data_transformed <- all_interval_data
all_interval_data_transformed[,c(6,13:83)] <- log(all_interval_data_transformed[,c(6,13:83)] + 1)


# -------------------------------------------------
# compare transformed and not-transformed data
# -------------------------------------------------

# Plot data
my_data <- all_interval_data[,c(1,2:17,32, 70:71, 75:81)]
my_data$rr_id <- gsub("_rr.*", "", my_data$rr_id)
a <- my_data %>%
  tidyr::pivot_longer(!rr_id) %>% 
  ggplot(aes(x=value, fill = rr_id, color = rr_id)) +
  facet_wrap(~ name, scales = "free") +
  geom_histogram(position="identity", alpha=0.4, bins = 50) +
  theme_classic()


# Splitting intergenic and genic data and coloring them differently
lala <- all_interval_data[,c(1,2:17,32, 70:71, 75:81)]
lala[, c(2:27)] <- log(lala[, c(2:27)] + 1)
lala$rr_id <- gsub("_rr.*", "", lala$rr_id)
b <- lala %>%
  pivot_longer(!rr_id) %>% 
  ggplot(aes(x=value, fill = rr_id, color = rr_id)) +
  facet_wrap(~ name, scales = "free") +
  geom_histogram(position="identity", alpha=0.4, bins = 50) +
  theme_classic()

# It looks like log transformations work nice
a/b

# filename = "~/git_projects/haplotype_v_snp_gwas/images/distributions_fixed_random_response.pdf",
ggsave(filename = "C:/Users/merri/git_projects/haplotype_v_snp_gwas/images/distributions_fixed_random_response.pdf",
       plot = a/b,
       width = 15,
       height = 16,
       units = "in")




# -----------------------------------------------------------------
#                     Plot generic things
# -----------------------------------------------------------------

# colnames(entropy_gwasCounts) <- gsub("shannon_entropy_", "", colnames(entropy_gwasCounts))
# colnames(entropy_gwasCounts) <- gsub("nam_", "NAM ", colnames(entropy_gwasCounts))
# colnames(entropy_gwasCounts) <- gsub("goodman_", "Goodman ", colnames(entropy_gwasCounts))
# colnames(entropy_gwasCounts) <- gsub("unfiltered_all", "Unfiltered Physiological", colnames(entropy_gwasCounts))
# colnames(entropy_gwasCounts) <- gsub("unfiltered_metabolite", "Unfiltered Metabolite", colnames(entropy_gwasCounts))
# colnames(entropy_gwasCounts) <- gsub("unfiltered_expression", "Unfiltered Expression", colnames(entropy_gwasCounts))
# colnames(entropy_gwasCounts) <- gsub("unfiltered_physiological", "Unfiltered Physiological", colnames(entropy_gwasCounts))
# colnames(entropy_gwasCounts) <- gsub("filtered_all", "Filtered Physiological", colnames(entropy_gwasCounts))
# colnames(entropy_gwasCounts) <- gsub("filtered_physiological", "Filtered Physiological", colnames(entropy_gwasCounts))
# colnames(entropy_gwasCounts) <- gsub("filtered_metabolite", "Filtered Metabolite", colnames(entropy_gwasCounts))
# colnames(entropy_gwasCounts) <- gsub("filtered_expression", "Filtered Expression", colnames(entropy_gwasCounts))


# -------------------------------------------
# historgrams of entropy/diversity by trait type
# -------------------------------------------

# set axis text size
axis_text_size <- 8

# Physiological NAM
temp <- reshape2::melt(entropy_gwasCounts[,1:3], id = c("Interval"))
pn <- ggplot(temp, aes(x = value, color = variable, fill = variable)) +
  geom_histogram(position="identity", alpha=0.4, bins = 50) +
  xlab("Shannon Entropy - NAM Physiological") +
  ylab("Count of Intervals") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), axis.title = element_text(size=axis_text_size), 
        legend.position="bottom", legend.title = element_text(size = 5), legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, 'cm')) +
  ylim(0, 20000)

# Physiological goodman
temp <- reshape2::melt(entropy_gwasCounts[,c(1,6,7)], id = c("Interval"))
pg <- ggplot(temp, aes(x = value, color = variable, fill = variable)) +
  geom_histogram(position="identity", alpha=0.4, bins = 50) +
  xlab("Shannon Entropy - Goodman Physiological") +
  ylab("Count of Intervals") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), axis.title = element_text(size=axis_text_size), 
        legend.position="bottom", legend.title = element_text(size = 5), legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, 'cm')) +
  ylim(0, 20000)


# Metabolite
temp <- reshape2::melt(entropy_gwasCounts[,c(1,4,5)], id = c("Interval"))
m <- ggplot(temp, aes(x = value, color = variable, fill = variable)) +
  geom_histogram(position="identity", alpha=0.4, bins = 50) +
  xlab("Shannon Entropy - Goodman Metabolite") +
  ylab("Count of Intervals") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), axis.title = element_text(size=axis_text_size), 
        legend.position="bottom", legend.title = element_text(size = 5), legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, 'cm')) +
  ylim(0, 20000)


# Expression
temp <- reshape2::melt(entropy_gwasCounts[,c(1,8,9)], id = c("Interval"))
e <- ggplot(temp, aes(x = value, color = variable, fill = variable)) +
  geom_histogram(position="identity", alpha=0.4, bins = 50) +
  xlab("Shannon Entropy - Goodman Expression") +
  ylab("Count of Intervals") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), axis.title = element_text(size=axis_text_size), 
        legend.position="bottom", legend.title = element_text(size = 5), legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, 'cm')) +
  ylim(0, 20000)

# Look at image together
# (pn + pg)/(m+e)

# Save to file
ggsave(
  filename = "/home/mbb262/git_projects/haplotype_v_snp_gwas/images/diversity_phys_metabolite.png",
  plot = (pn + pg)/(m+e),
  width = 6.5,
  height = 3.5,
  units = "in",
  dpi = 300
)


# -----------------------------------------------------------------
#                         Run models
# -----------------------------------------------------------------

# RNA models
# entropy ~ RNA expression + #gerp + range_type +(1|LD) + (1|#counts) + (1|interval size)

# different sets of models to run
# nam filtered/unfiltered physiological - all (np_all)
# nam filtered/unfiltered physiological - genic (np_genic)
# nam filtered/unfiltered physiological - intergenic (np_intergenic)

# goodman filtered/unfiltered physiological - all (gfp_all)
# goodman filtered/unfiltered physiological - genic (gfp_genic)
# goodman filtered/unfiltered physiological - intergenic (gfp_intergenic)

# goodman filtered/unfiltered metabolite - all (gm_all)
# goodman filtered/unfiltered metabolite - genic (gm_genic)
# goodman filtered/unfiltered metabolite - intergenic (gm_intergenic)

# goodman filtered/unfiltered expression - all (ge_all)
# goodman filtered/unfiltered expression - genic (ge_genic)
# goodman filtered/unfiltered expression - intergenic (ge_intergenic)



# -------------------------
#   No model, just plot
# entropy vs expression
# -------------------------

# scatter plot
ggplot(all_interval_data_transformed, aes(x = Mature_Leaf_8_rna, y = shannon_entropy_nam_filtered_all, color = range_type)) + 
  geom_smooth(method = "lm") +
  geom_point(alpha = 0.5)+
  labs(y="Shannon Entropy - Physiological Traits", x = "Walley Mature Leaf Expression") +
  geom_text(x=3, y=30, label="Scatter plot")

# hex2bin
all_interval_data_transformed %>% 
  dplyr::filter(Mature_Leaf_8_rna & shannon_entropy_nam_filtered_all) %>% 
  ggplot() +
  aes(x = Mature_Leaf_8_rna, y = shannon_entropy_nam_filtered_all) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis", name = "Count") +
  xlab("Walley Max Expression Across Tissues") 


# -----------------------------------------------------------------------
#                             NAM Models
# -----------------------------------------------------------------------

# --------------------
# rearrange data
# --------------------

np_all_data <- tidyr::pivot_longer(all_interval_data_transformed[,c(1,2,6,13, 15,16,23,24,65,82)], 
                                   cols = shannon_entropy_nam_filtered_all:shannon_entropy_nam_unfiltered_all,
                                   names_to = "entropy_filtering_level",
                                   values_to = "entropy")
np_all_data <- tidyr::pivot_longer(np_all_data,
                                   cols = sum_gwa_hitsnam_filtered_all:sum_gwa_hitsnam_unfiltered_all,
                                   names_to = "sum_filtering_level",
                                   values_to = "sum_gwa_hits")

# Clean up filtered and unfiltered names for plotting later on
np_all_data$entropy_filtering_level <- gsub("shannon_entropy_nam_f", "F", np_all_data$entropy_filtering_level)
np_all_data$entropy_filtering_level <- gsub("shannon_entropy_nam_u", "U", np_all_data$entropy_filtering_level)
np_all_data$entropy_filtering_level <- gsub("_all", "", np_all_data$entropy_filtering_level)

# generic plot
ggplot(np_all_data, aes(x = Mature_Leaf_8_rna, y = entropy, color = entropy_filtering_level)) + 
  geom_smooth(method = "lm") +
  geom_point(alpha = 0.5)+
  labs(y="Shannon Entropy - Physiological Traits", x = "Walley Mature Leaf Expression") +
  geom_text(x=3, y=30, label="Scatter plot")


# --------------------
# nam filtered/unfiltered physiological all fixed effects - all (np_all_fixed)
# --------------------

# run model - all fixed effects
np_all_fixed <- lm(entropy ~ Mature_Leaf_8_rna  + 
                     range_type + 
                     entropy_filtering_level+ 
                     interval_size + 
                     sum_gwa_hits + 
                     average_r2_ld_nam, 
                   data = np_all_data)
summary(np_all_fixed) # model is singular
anova(np_all_fixed)
np_all_data$predicted <- predict(np_all_fixed, re.form=NA) # predict values
np_all_fixed_plot <- ggplot(np_all_data, aes(x = Mature_Leaf_8_rna, y = predicted, color = entropy_filtering_level))+
  geom_smooth(method = "lm") +
  geom_point(alpha = 0.5) +
  labs(y="Fitted Shannon Entropy - NAM physiological", 
       x = "Walley Mature Leaf Expression", 
       title = "NAM Physiological, all fixed effects")


# --------------------
# nam filtered/unfiltered physiological - all (np_all)
# --------------------

# run model
np_all <- lme4::lmer(entropy ~ Mature_Leaf_8_rna  + 
                       range_type + filtering_level+ (1|interval_size) + (1|sum_gwa_hitnam_filtered_all) + 
                       (1|average_r2_ld_nam), 
                     data = np_all_data)
isSingular(np_all)
summary(np_all) # model is singular
anova(np_all)
np_all_data$predicted <- predict(np_all, re.form=NA) # predict values
ggplot(np_all_data, aes(x = Mature_Leaf_8_rna, y = predicted))+
  geom_point()


# --------------------
# nam filtered/unfiltered physiological - genic (np_genic)
# --------------------

np_genic_data <- np_all_data %>%  filter(range_type == "genic")
np_genic <- lm(entropy ~ Mature_Leaf_8_rna + gerp_count_thresh_level_6 + 
                         interval_size + sum_gwa_hits + 
                         average_r2_ld_nam, 
                       data = np_genic_data)
isSingular(np_genic)
summary(np_genic) # model is singular
anova(np_genic)
np_genic_data$predicted <- predict(np_genic, re.form=NA)
ggplot(np_genic_data, aes(x = Mature_Leaf_8_rna, y = predicted, color = entropy_filtering_level))+
  geom_smooth(method = "lm") +
  geom_point(alpha = 0.5) +
  labs(y="Fitted Shannon Entropy - NAM physiological", 
       x = "Walley Mature Leaf Expression", 
       title = "NAM Physiological - Genic Regions, all fixed effects")



# --------------------
# nam filtered/unfiltered physiological - intergenic (np_intergenic)
# --------------------


# -----------------------------------------------------------------------
#                  Goodman physiological models
# -----------------------------------------------------------------------

gp_all_data <- tidyr::pivot_longer(all_interval_data_transformed[,c(1,2,6,14,19,20,27,28,65,82)], 
                                   cols = shannon_entropy_goodman_filtered_physiological:shannon_entropy_goodman_unfiltered_physiological,
                                   names_to = "entropy_filtering_level",
                                   values_to = "entropy")
gp_all_data <- tidyr::pivot_longer(gp_all_data,
                                   cols = sum_gwa_hitsgoodman_filtered_physiological:sum_gwa_hitsgoodman_unfiltered_physiological,
                                   names_to = "sum_filtering_level",
                                   values_to = "sum_gwa_hits")

# generic plot
ggplot(gp_all_data, aes(x = Mature_Leaf_8_rna, y = entropy, color = entropy_filtering_level)) + 
  geom_smooth(method = "lm") +
  geom_point(alpha = 0.5)+
  labs(y="Shannon Entropy - Physiological Traits", x = "Walley Mature Leaf Expression") +
  geom_text(x=3, y=30, label="Scatter plot")


# --------------------
# goodman filtered/unfiltered physiological all fixed effects - all (gp_all_fixed)
# --------------------

# run model - all fixed effects
gp_all_fixed <- lm(entropy ~ Mature_Leaf_8_rna  + 
                     range_type + 
                     entropy_filtering_level+ 
                     interval_size + 
                     sum_gwa_hits + 
                     average_r2_ld_282, 
                   data = gp_all_data)
summary(gp_all_fixed) # model is singular
anova(gp_all_fixed)
gp_all_data$predicted <- predict(gp_all_fixed, re.form=NA) # predict values
gp_all_fixed_plot <- ggplot(gp_all_data, aes(x = Mature_Leaf_8_rna, y = predicted, color = entropy_filtering_level))+
  geom_smooth(method = "lm") +
  geom_point(alpha = 0.5) +
  labs(y="Fitted Shannon Entropy - goodman physiological", 
       x = "Walley Mature Leaf Expression", 
       title = "goodman Physiological, all fixed effects")


# -----------------------------------------------------------------------
#                  Goodman metabolite models
# -----------------------------------------------------------------------


# -----------------------------------------------------------------------
#                  Goodman expression models
# -----------------------------------------------------------------------

#melt dataframe twice to get the approporate structue
ge_all_data <- tidyr::pivot_longer(all_interval_data_transformed[,c(1,2,6,14, 21,22,29,30,65,82)], 
                                   cols = shannon_entropy_goodman_filtered_expression:shannon_entropy_goodman_unfiltered_expression,
                                   names_to = "entropy_filtering_level",
                                   values_to = "entropy")
ge_all_data <- tidyr::pivot_longer(ge_all_data,
                                   cols = sum_gwa_hitsgoodman_filtered_expression:sum_gwa_hitsgoodman_unfiltered_expression,
                                   names_to = "sum_filtering_level",
                                   values_to = "sum_gwa_hits")

# generic plot
ggplot(ge_all_data, aes(x = Mature_Leaf_8_rna, y = entropy, color = entropy_filtering_level)) + 
  geom_smooth(method = "lm") +
  geom_point(alpha = 0.5)+
  labs(y="Shannon Entropy - Physiological Traits", x = "Walley Mature Leaf Expression") +
  geom_text(x=3, y=30, label="Scatter plot")


# --------------------
# goodman filtered/unfiltered expression - only fixed effects - all (ge_all_fixed) 
# --------------------

# run model - all fixed effects
ge_all_fixed <- lm(entropy ~ Mature_Leaf_8_rna  + 
                     range_type + 
                     entropy_filtering_level+ 
                     interval_size + 
                     sum_gwa_hits + 
                     average_r2_ld_282, 
                   data = ge_all_data)
summary(ge_all_fixed) # model is singular
anova(ge_all_fixed)
ge_all_data$predicted <- predict(ge_all_fixed, re.form=NA) # predict values
ge_all_fixed_plot <- ggplot(ge_all_data, aes(x = Mature_Leaf_8_rna, y = predicted, color = entropy_filtering_level))+
  geom_smooth(method = "lm") +
  geom_point() +
  labs(y="Fitted Shannon Entropy - 282 Expression", 
       x = "Walley Mature Leaf Expression", 
       title = "Goodman Expression, all fixed effects")

ggsave(ge_all_fixed_plot,"/workdir/mbb262/ge_all_fixed.png")


# --------------------
# goodman filtered/unfiltered expression - all (ge_all)
# --------------------

# run model - fixed and random
ge_all <- lme4::lmer(entropy ~ Mature_Leaf_8_rna  + 
                       range_type + 
                       entropy_filtering_level+ 
                       (1|interval_size) + 
                       (1|sum_gwa_hits) + 
                       (1|average_r2_ld_282), 
                     data = ge_all_data)
isSingular(ge_all)
summary(ge_all) # model is singular
anova(ge_all)
ge_all_data$predicted <- predict(ge_all, re.form=NA) # predict values
ge_all_plot <- ggplot(ge_all_data, aes(x = Mature_Leaf_8_rna, y = predicted, color = entropy_filtering_level))+
  geom_smooth(method = "lm") +
  geom_point() +
  labs(y="Fitted Shannon Entropy - 282 Expression", 
       x = "Walley Mature Leaf Expression",
       title = "Goodman Expression, fixed and random effects")


# --------------------
# goodman filtered/unfiltered expression - genic (ge_genic)
# --------------------

ge_genic_data <- ge_all_data %>% filter(range_type == "genic")

# run model - fixed and random
ge_genic <- lme4::lmer(entropy ~ Mature_Leaf_8_rna  + 
                         entropy_filtering_level+ 
                         (1|interval_size) + 
                         (1|sum_gwa_hits) + 
                         (1|average_r2_ld_282), 
                       data = ge_genic_data)
isSingular(ge_genic)
summary(ge_genic) # model is singular
anova(ge_genic)
ge_genic_data$predicted <- predict(ge_genic, re.form=NA) # predict values
ge_genic_plot <- ggplot(ge_genic_data, aes(x = Mature_Leaf_8_rna, y = predicted))+
  geom_point()

# --------------------
# goodman filtered/unfiltered expression - intergenic (ge_intergenic)
# --------------------



