# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2024-01-19
# Updated... 2024-01-19
#
# Description:
# Model TE families in individual models JUST FOR B73 to see if results are
# consistent between the inbreds and hybrids
#
#  Expression ij ~ Total TE bp (within X kb) at 12 distances ij + 
#                  12 UMRs + mean PhyloP + mean dnds + 
#                  mean B73 expression + tissue
# (Where i is all genes and j is each TE)
# ------------------------------------------------------------------------------

# Load packages
library(dplyr)
library(ggplot2)
library(foreach)
library(doParallel)
library(data.table)

# Set global variables
num_cores <- 85


# Data formatting --------------------------------------------------------------

# Load data (no phenotype version)
all_model_data_noTEs <- data.table::fread("/workdir/mbb262/model_data/all_merged_model_data_noTEs.csv", nThread = 40)

# Only use complete cases
all_model_data_noTEs <- all_model_data_noTEs %>% 
  filter(complete.cases(.) & age != "58 days" & genotype == "B73")

# Load in TE data
te <- data.table::fread("/workdir/mbb262/teFams_liftoffgene_distbins.2023-12-22.txt", nThread = 40)

# Rearrange columns and rename
te <- te %>% select(genotype, X, fam, F_500, F_1000, F_1500, F_3000, F_4500, F_6000, F_7500, F_9000, F_10500, F_12000, F_13500, F_15000)
te$genotype <- toupper(te$genotype)
colnames(te) <- c("genotype", "gene", "family", "tebp_1_500", "tebp_501_1000", "tebp_1001_1500", 
                  "tebp_1501_3000", "tebp_3001_4500", "tebp_4501_6000", "tebp_6001_7500", 
                  "tebp_7501_9000", "tebp_9001_10500", "tebp_10501_12000", 
                  "tebp_12001_13500", "tebp_13501_15000")
te[1:5,1:5]

# Threshold values according to window size
# TEMPORARY 
# Threshold TE counts, counts larger than a window should just = window
te$tebp_1_500[te$tebp_1_500 > 500] <- 500
te$tebp_501_1000[te$tebp_501_1000 > 500] <- 500
te$tebp_1001_1500[te$tebp_1001_1500 > 500] <- 500
te$tebp_1501_3000[te$tebp_1501_3000 > 1500] <- 1500
te$tebp_3001_4500[te$tebp_3001_4500 > 1500] <- 1500
te$tebp_4501_6000[te$tebp_4501_6000 > 1500] <- 1500
te$tebp_6001_7500[te$tebp_6001_7500 > 1500] <- 1500
te$tebp_7501_9000[te$tebp_7501_9000 > 1500] <- 1500
te$tebp_9001_10500[te$tebp_9001_10500 > 1500] <- 1500
te$tebp_10501_12000[te$tebp_10501_12000 > 1500] <- 1500
te$tebp_12001_13500[te$tebp_12001_13500 > 1500] <- 1500
te$tebp_13501_15000[te$tebp_13501_15000 > 1500] <- 1500

# Scale by window size
te$scaledTEBP_1_500 <- te$tebp_1_500/500
te$scaledTEBP_501_1000 <- te$tebp_501_1000/500
te$scaledTEBP_1001_1500 <- te$tebp_1001_1500/500
te$scaledTEBP_1501_3000 <- te$tebp_1501_3000/1500
te$scaledTEBP_3001_4500 <- te$tebp_3001_4500/1500
te$scaledTEBP_4501_6000  <- te$tebp_4501_6000/1500
te$scaledTEBP_6001_7500 <- te$tebp_6001_7500/1500
te$scaledTEBP_7501_9000 <- te$tebp_7501_9000/1500
te$scaledTEBP_9001_10500 <- te$tebp_9001_10500/1500
te$scaledTEBP_10501_12000 <- te$tebp_10501_12000/1500
te$scaledTEBP_12001_13500 <- te$tebp_12001_13500/1500
te$scaledTEBP_13501_15000 <- te$tebp_13501_15000/1500

# Subset only the columns needed
te <- te[,c("genotype", "gene", "family", "scaledTEBP_1_500", "scaledTEBP_501_1000",
            "scaledTEBP_1001_1500", "scaledTEBP_1501_3000", "scaledTEBP_3001_4500",
            "scaledTEBP_4501_6000", "scaledTEBP_6001_7500", "scaledTEBP_7501_9000",
            "scaledTEBP_9001_10500", "scaledTEBP_10501_12000", "scaledTEBP_12001_13500",
            "scaledTEBP_13501_15000")]


# Helper functions -------------------------------------------------------------

modelDistance <- function(variableName, simplifiedName, teFamilyName, df) {
  # Create model statement
  modelStatement <- as.formula(paste0("expression ~ ", variableName, " + mean_expression +  mean_phylop + mean_dnds + age + tissue"))
  
  # Run model
  runModel <- lm(modelStatement, data = df)
  
  # Collect model metrics
  individualFullModelStats <- data.frame(adj_r2 = broom::glance(runModel)$adj.r.squared,
                                         pvalue = broom::glance(runModel)$p.value,
                                         model = paste0(simplifiedName, " Model"),
                                         family = teFamilyName)
  
  # Collect term effects
  fullModelnoUMREffects <- broom::tidy(runModel)
  fullModelnoUMREffects$model <- rep(paste0(simplifiedName, " Model"), rep(nrow(fullModelnoUMREffects)))
  fullModelnoUMREffects$family <- rep(teFamilyName, rep(nrow(fullModelnoUMREffects)))
  
  # Save both as list
  return(list(individualFullModelStats, fullModelnoUMREffects))
}


# Run models -------------------------------------------------------------------

# Find things to iterate over
uniqueTEfamilies <- unique(te$family)
length(uniqueTEfamilies)

# Register parallel backend
registerDoParallel(cores = num_cores)

# Loop through each family individually to 
results_list <- foreach(i = 1:length(uniqueTEfamilies), .packages = c("data.table","broom")) %dopar% {
  
  # Subset te matrix to just that family
  sub_te <- te[family == uniqueTEfamilies[i]]
  
  # Merge with phenotypic data (keyed merging saved on .5 sec per te family)
  setkeyv(all_model_data_noTEs,  c("genotype", "gene"))
  setkeyv(sub_te,  c("genotype", "gene"))
  sub_all_model_data <- merge(all_model_data_noTEs, sub_te, all.x = TRUE)
  
  # There are some rare instances where a TE family is only present next to one
  # gene. If this gene is not expressed, it will not be present in sub_all_model_data
  # And merging will result in a dataframe with nothing in it. These types of TEs 
  # are discard from the model later. For now, an if statemement is below to yeet 
  # out of the foreach statement
  if(length(unique(sub_all_model_data$family)) > 1){
    # Repeat the missing TE family values
    sub_all_model_data$family <- rep(uniqueTEfamilies[i], nrow(sub_all_model_data))
    
    # Replace NAs with zeros in the distances
    sub_all_model_data[, 23:34][is.na(sub_all_model_data[, 23:34])] <- 0
    
    # Run model
    model1 <- lm(expression ~ scaledTEBP_1_500 + scaledTEBP_501_1000 + scaledTEBP_1001_1500 + scaledTEBP_1501_3000 + scaledTEBP_3001_4500 + scaledTEBP_4501_6000 + scaledTEBP_6001_7500 + scaledTEBP_7501_9000 + scaledTEBP_9001_10500 + scaledTEBP_10501_12000 + scaledTEBP_12001_13500 + scaledTEBP_13501_15000 + mean_expression + mean_phylop + mean_dnds + age + tissue,
                 data = sub_all_model_data)
    
    # Put model summary in a cleaned up df, add metadata
    jointModelnoUMREffects <- broom::tidy(model1)
    jointModelnoUMREffects$model <- rep("Joint Term Model", rep(nrow(jointModelnoUMREffects)))
    jointModelnoUMREffects$family <- rep(uniqueTEfamilies[i], rep(nrow(jointModelnoUMREffects)))
    
    # Collect overall model statistics, add metadata
    jointModelStats <- data.frame(adj_r2 = broom::glance(model1)$adj.r.squared, 
                                  pvalue = broom::glance(model1)$p.value, 
                                  model = "Joint Term Model",
                                  family = uniqueTEfamilies[i])
    
    # For each individual family, use outside function
    m1 <- modelDistance("scaledTEBP_1_500", "1-500", uniqueTEfamilies[i], sub_all_model_data)
    m2 <- modelDistance("scaledTEBP_501_1000", "501_1000", uniqueTEfamilies[i], sub_all_model_data)
    m3 <- modelDistance("scaledTEBP_1001_1500", "1001_1500", uniqueTEfamilies[i], sub_all_model_data)
    m4 <- modelDistance("scaledTEBP_1501_3000", "1501_3000", uniqueTEfamilies[i], sub_all_model_data)
    m5 <- modelDistance("scaledTEBP_3001_4500", "3001_4500", uniqueTEfamilies[i], sub_all_model_data)
    m6 <- modelDistance("scaledTEBP_4501_6000", "4501_6000", uniqueTEfamilies[i], sub_all_model_data)
    m7 <- modelDistance("scaledTEBP_6001_7500", "6001_7500", uniqueTEfamilies[i], sub_all_model_data)
    m8 <- modelDistance("scaledTEBP_7501_9000", "7501_9000", uniqueTEfamilies[i], sub_all_model_data)
    m9 <- modelDistance("scaledTEBP_9001_10500", "9001_10500", uniqueTEfamilies[i], sub_all_model_data)
    m10 <- modelDistance("scaledTEBP_10501_12000", "10501_12000", uniqueTEfamilies[i], sub_all_model_data)
    m11 <- modelDistance("scaledTEBP_12001_13500", "12001_13500", uniqueTEfamilies[i], sub_all_model_data)
    m12 <- modelDistance("scaledTEBP_13501_15000", "13501_15000", uniqueTEfamilies[i], sub_all_model_data)
    
    # Combine partial effect estimates
    all_tefamily_model_effects <- rbind(jointModelnoUMREffects,
                                        m1[[2]], m2[[2]], m3[[2]], m4[[2]], m5[[2]], m6[[2]], 
                                        m7[[2]], m8[[2]], m9[[2]], m10[[2]], m11[[2]], m12[[2]])
    
    # Full model stats
    allModelSummaries <- rbind(jointModelStats, 
                               m1[[1]], m2[[1]], m3[[1]], m4[[1]], m5[[1]], m6[[1]], 
                               m7[[1]], m8[[1]], m9[[1]], m10[[1]], m11[[1]], m12[[1]])
    
    # Return all results
    list(all_tefamily_model_effects, allModelSummaries)
  }
}  


# Format results ---------------------------------------------------------------

# Check lengths of outputs
length(uniqueTEfamilies)
length(results_list)

# Split into separate dfs
modelEffectResults <- do.call(rbind, lapply(results_list, `[[`, 1))
modelSummaryResults <- do.call(rbind, lapply(results_list, `[[`, 2))

# Write to file
data.table::fwrite(modelEffectResults, "/workdir/mbb262/model_data/model_results/inbred_individualFamilies/inbred_B73_individual_te_family_model_term_effects_by_distance.csv", nThread = 40)
data.table::fwrite(modelSummaryResults, "/workdir/mbb262/model_data/model_results/inbred_individualFamilies/inbred_B73_individual_te_family_model_wide_effects_by_distance.csv", nThread = 40)


# Filter families and plot results ---------------------------------------------


# Load packages
library(dplyr)
library(ggplot2)
library(data.table)


# Load in TE info --------------------------------------------------------------

# Merge TE family names
teFamilyNames <- read.delim("/workdir/mbb262/tefamname_crossreference_for_merritt.txt")
class2Super <- read.delim("/workdir/mbb262/TE_content_across_NAM_genotypes_by_collapsedfam.2023-07-14.txt")
teFamilyNames <- merge(teFamilyNames, class2Super, by = "collapsedFam")

# Format TE counts to filter by
te <- data.table::fread("/workdir/mbb262/teFams_liftoffgene_distbins.2023-12-22.txt", nThread = 40)

# Rearrange columns and rename
te <- te %>% select(genotype, X, fam, F_500, F_1000, F_1500, F_3000, F_4500, F_6000, F_7500, F_9000, F_10500, F_12000, F_13500, F_15000)
te$genotype <- toupper(te$genotype)
colnames(te) <- c("genotype", "gene", "te_family", "tebp_1_500", "tebp_501_1000", "tebp_1001_1500", 
                  "tebp_1501_3000", "tebp_3001_4500", "tebp_4501_6000", "tebp_6001_7500", 
                  "tebp_7501_9000", "tebp_9001_10500", "tebp_10501_12000", 
                  "tebp_12001_13500", "tebp_13501_15000")

# Threshold values according to window size --> TEMPORARY 
# Threshold TE counts, counts larger than a window should just = window
te$tebp_1_500[te$tebp_1_500 > 500] <- 500
te$tebp_501_1000[te$tebp_501_1000 > 500] <- 500
te$tebp_1001_1500[te$tebp_1001_1500 > 500] <- 500
te$tebp_1501_3000[te$tebp_1501_3000 > 1500] <- 1500
te$tebp_3001_4500[te$tebp_3001_4500 > 1500] <- 1500
te$tebp_4501_6000[te$tebp_4501_6000 > 1500] <- 1500
te$tebp_6001_7500[te$tebp_6001_7500 > 1500] <- 1500
te$tebp_7501_9000[te$tebp_7501_9000 > 1500] <- 1500
te$tebp_9001_10500[te$tebp_9001_10500 > 1500] <- 1500
te$tebp_10501_12000[te$tebp_10501_12000 > 1500] <- 1500
te$tebp_12001_13500[te$tebp_12001_13500 > 1500] <- 1500
te$tebp_13501_15000[te$tebp_13501_15000 > 1500] <- 1500
te[1:5,1:5]

# Filter to only TE families present in at least 10 genes (and not alleles)
countFamilies <- te %>%
  select(gene, te_family) %>%
  distinct() %>%
  group_by(te_family) %>%
  summarise(n = n()) %>%
  filter(n > 10)
nrow(countFamilies) #after filter
length(unique(te$te_family)) # before filter

# Collect families present in at least 25 taxa-gene-distance combos
# (with at least 1 base pair in each distance category)
# For the joint term models, all TE distances must pass filters
meltTE <- reshape2::melt(te, id.vars = c("genotype", "gene", "te_family"))
keepableTEDistances_single <- meltTE %>%
  filter(value > 0) %>% # Filter out distances with no te base pairs
  group_by(te_family, variable) %>%
  summarise(n = n()) %>%
  group_by(te_family) %>%
  filter(all(n > 25)) %>% # Filter out TE families that are in each distance less than 25 times - this value removes some of the major outliers
  filter(te_family != "FAM14838") # Filter out one problem family with extreme effects that passes every filter in the inbreds
keepableTEDistances_single$dummyFilterVar <- paste0(keepableTEDistances_single$te_family, "_", keepableTEDistances_single$variable)
nrow(keepableTEDistances_single)
length(unique(te$te_family))


# Format Inbred Model Effects --------------------------------------------------

# Load in full model stats, Merge model results with TE family names, filter
inbred_fullModelResults <- data.table::fread("/workdir/mbb262/model_data/model_results/inbred_individualFamilies/inbred_B73_individual_te_family_model_wide_effects_by_distance.csv")
inbred_fullModelResults <- merge(inbred_fullModelResults, teFamilyNames, by.x = "family", by.y = "dumbindex") %>% 
  select(-pvalue) 
# %>% 
#   filter(model != "Joint Term Model")

# Load in model effects, turn effect values numeric, create a summary variable to filter out distances
inbred_allModelEffects <- data.table::fread("/workdir/mbb262/model_data/model_results/inbred_individualFamilies/inbred_B73_individual_te_family_model_term_effects_by_distance.csv")

# Remove model levels
removeLevels <- c("(Intercept)", "mean_expression", "mean_phylop", "mean_dnds", "age42 days", "age50 days", "age58 days", "tissueMiddle section of adult flag leaf", "tissueV3 leaf base", "tissueV3 leaf tip")
inbred_allModelEffects <- inbred_allModelEffects %>% 
  filter(!term %in% removeLevels) %>% 
  # filter(model != "Joint Term Model") %>% 
  select(-std.error, -statistic)

# Merge both together
inbred <- merge(inbred_fullModelResults, inbred_allModelEffects, by = c("family", "model"))

# Make data for plotting
inbred$Term <- gsub("scaledTEBP_", "", inbred$term)
inbred$Term <- gsub("_", "-", inbred$Term)

# Create a dummy variable to help with filtering
inbred$dummyFilterVar <- paste0(inbred$family, "_", gsub("scaledTEBP_", "tebp_", inbred$term))

# Subset results based on the two vectors/criteria then combine
inbred <- inbred %>%
  filter(dummyFilterVar %in% keepableTEDistances_single$dummyFilterVar) %>%
  filter(family %in% countFamilies$te_family)

# Number of unique inbreds
length(unique(inbred$family)) # 1605

# Change ordering of factors
inbred$Term <- as.character(inbred$Term)
inbred$Term <- factor(inbred$Term, 
                      levels=c("1-500","501-1000", "1001-1500", 
                               "1501-3000", "3001-4500", "4501-6000",
                               "6001-7500", "7501-9000", "9001-10500", 
                               "10501-12000", "12001-13500", "13501-15000"))


# Plot inbred Data -------------------------------------------------------------

# Plot R2 against TE family copy number ---------------------------

# COmbine with superfamily information
supFam <- read.delim("/workdir/mbb262/TE_content_across_PHZ51_by_collapsedfam.2023-07-14.txt") %>% select(-bp, -Classification)

joint_inbred_fullModelResults <- inbred %>% filter(model == "Joint Term Model") %>% select(family, adj_r2, bp, collapsedFam) %>% distinct()
joint_inbred_fullModelResults <- merge(joint_inbred_fullModelResults, supFam, by = "collapsedFam")

# Make a threshold line
temp <- quantile(joint_inbred_fullModelResults$adj_r2, .99)
temp

# Plot across super-families
r2ByTotalBp <- ggplot(joint_inbred_fullModelResults, aes(x = log10(bp), y = adj_r2, color = sup)) +
  geom_point() +
  xlab("log10(Total TE Family Base Pairs)") +
  ylab("Adjusted R2") +
  theme_classic() +
  geom_hline(yintercept = temp) +
  theme(axis.text=element_text(size=axis_text_size),
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size),
        legend.position="bottom",
        legend.text = element_text(size = axis_text_size),
        legend.title=element_text(size = legend_text_size),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) 
r2ByTotalBp

# Get the top 5 families
joint_inbred_fullModelResults %>%
  slice_max(adj_r2, n = 10)

ggsave("/home/mbb262/git_projects/te_ase_nam/figs/joint_inbred_B73Only_modelR2_filteredByTEPrevalance_totalTEBP.png",
       plot = r2ByTotalBp,
       width = 6.5,
       height = 5.5,
       units = "in",
       dpi = "retina")


# Plot what michelle suggested -------------------------------------------------

# What michelle suggested:
# 1) grouping by family size (<1Mb, 1-10Mb, 10-100Mb, >100Mb) by effect

# 2) proportion of the family inserted in these 15kb near genes (bp summed in 15kb window / total family bp) by effect
# might help test family trends? Would a lineplot (with heavy alpha to overlap) be 
# helpful to show if the outliers are consistent outliers across the distance range?

#1 --------------------------------
# Cut values into bins
inbred$categories <- cut(inbred$bp, breaks = c(0, 1000000, 10000000, 100000000, Inf),
                                             labels = c("<1Mb", "1-10Mb", "10-100Mb", ">100Mb"))

# Subset to only joint effect models
sub_inbred_allModelEffects_joint <- inbred %>% 
  filter(model == "Joint Term Model" )

# Plot against effect
famSize <- ggplot(sub_inbred_allModelEffects_joint, aes(y = estimate, x = Term, fill = categories)) +
  geom_boxplot() +
  xlab("Distance from Gene") +
  ylab("Effect") +
  labs(fill = "TE Family Size") +
  ggtitle("B73 only model with coord_cartesian(ylim = c(-15,15)")+
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette)+
  coord_cartesian(ylim = c(-15,15), expand = F)
famSize

ggsave("/home/mbb262/git_projects/te_ase_nam/figs/joint_inbred_B73Only_model_teFamilySize_by_effect_boxplot.png",
       plot = famSize,
       width = 6.5,
       height = 5.5,
       units = "in",
       dpi = "retina")


#2 --------------------------------

# Calculate the proportion of the family inserted in these 15kb near genes (bp summed in 15kb window / total family bp)

# Calculate sum of TE bp within 15 kb of a gene for each te family
bp_summed <- meltTE %>% 
  group_by(te_family) %>% 
  summarise(sum = sum(value))

# Merge together
together <- merge(bp_summed, sub_inbred_allModelEffects_joint, by.x  = "te_family", by.y = "family")

# Values over 1
together$propFamily15kb <- together$sum/together$bp
summary(together$propFamily15kb)

# Cut values into bins
together$propBins <- cut(together$propFamily15kb, breaks = c(0, 0.25, 0.50, 0.75, Inf),
                         labels = c("<25%", "25-50%", "50-75%", ">75%"))

proportion_family15_box <- ggplot(together, aes(x = Term, y= estimate, color = propBins)) +
  geom_boxplot() +
  xlab("Distance from Gene") +
  ylab("Effect") +
  labs(fill = "Proportion of TE Family within 15Kb") +
  ggtitle("B73 only model with coord_cartesian(ylim = c(-15,15)")+
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette)+
  coord_cartesian(ylim = c(-15,15), expand = F)
proportion_family15_box

ggsave("/home/mbb262/git_projects/te_ase_nam/figs/joint_inbred_B73Only_model_teFamilySize_by_effect_boxplot_binnedDistance.png",
       plot = proportion_family15_box,
       width = 6.5,
       height = 6.5,
       units = "in",
       dpi = "retina")


















