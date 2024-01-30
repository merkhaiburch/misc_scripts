# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-10-23
# Updated... 2023-12-21
#
# Description:
# Run a series of models from simple to complex looking at the relationship between
# gene expression and TE count or the number of TE base pairs
#
# Model 2
#  Expression ij ~ Total TE bp (within X kb) at 12 distances ij + 
#                  UMRs at 12 distances + mean PhyloP + mean dnds + 
#                  mean B73 expression + tissue
# (Where i is all genes and j is each TE)
#
# NEEDS WORK!!!!!!! I UPDATED THE SOURCE FILES AND THIGS NEED TO BE CHANGED
# ------------------------------------------------------------------------------

# Load packages
library(dplyr)
library(ggplot2)
library(ggpubr)

# Get data from blfs1
# scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/te_ase_nam/model_data/all_merged_model_data.csv ./

# Load data
all_model_data <- data.table::fread("/workdir/mbb262/model_data/all_merged_model_data.csv", nThread = 40)

# Only use complete cases
all_model_data <- all_model_data %>% 
  filter(complete.cases(.) & age != "58 days" & genotype == "B73")


## Joint/full term model -------------------------------------------------------
#  Expression ij ~ Total TE bp (within X kb) at 12 distances ij + 
#                  UMRs at 12 distances + mean PhyloP + mean dnds + 
#                  mean B73 expression + tissue
# (Where i is all genes and j is each TE)

# Scaled with UMRs
model2 <- lm(expression ~ scaledTEBP_1_500 + scaledTEBP_501_1000 + scaledTEBP_1001_1500 + scaledTEBP_1501_3000 + scaledTEBP_3001_4500 + scaledTEBP_4501_6000 + scaledTEBP_6001_7500 + scaledTEBP_7501_9000 + scaledTEBP_9001_10500 + scaledTEBP_10501_12000 + scaledTEBP_12001_13500 + scaledTEBP_13501_15000 + umrCount_1_500 + umrCount_1001_1500 + umrCount_10501_12000 + umrCount_12001_13500 + umrCount_13501_15000 + umrCount_1501_3000 + umrCount_3001_4500 + umrCount_4501_6000 + umrCount_501_1000 + umrCount_6001_7500 + umrCount_7501_9000 + umrCount_9001_10500 + mean_expression + mean_phylop + mean_dnds + tissue + age,
             data = all_model_data)
summary(model2)

# Put model summary in a cleaned up df, add metadata
jointModelnoUMREffects <- broom::tidy(model2)
jointModelnoUMREffects$model <- rep("Joint Term Model", rep(nrow(jointModelnoUMREffects)))

# Collect model statistics, add metadata
jointModelStats <- data.frame(adj_r2 = broom::glance(model2)$adj.r.squared, 
                              pvalue = broom::glance(model2)$p.value, 
                              model = "Joint Term Model")


# Calculate the percent variance explained by each term ------------------------
# R^2(full model) - R^2(model without the term)

# Create a df to collect the results
terms_cycle_through <- c("scaledTEBP_1_500", "scaledTEBP_501_1000", "scaledTEBP_1001_1500", 
                         "scaledTEBP_1501_3000", "scaledTEBP_3001_4500", "scaledTEBP_4501_6000", 
                         "scaledTEBP_6001_7500", "scaledTEBP_7501_9000", "scaledTEBP_9001_10500", 
                         "scaledTEBP_10501_12000", "scaledTEBP_12001_13500", "scaledTEBP_13501_15000", 
                         "umrCount_1_500", "umrCount_1001_1500", "umrCount_10501_12000", 
                         "umrCount_12001_13500", "umrCount_13501_15000", "umrCount_1501_3000", 
                         "umrCount_3001_4500", "umrCount_4501_6000", "umrCount_501_1000", 
                         "umrCount_6001_7500", "umrCount_7501_9000", "umrCount_9001_10500",
                         "mean_expression", "mean_phylop", "mean_dnds", "tissue", "age")

# Calculate the R2 value of each variable after dropping each term from the model
# Run full model, collect R2
fullRegressionFormula <- as.formula(paste("expression ~ ", paste(terms_cycle_through, collapse="+")))
full_model <- lm(formula = fullRegressionFormula, data = all_model_data)
full_model_r2 <- broom::glance(full_model)$adj.r.squared

# Outside loop result collector
outsideLoopTermResult <- c()

# Run individual models dropping a single term
# loop over terms
for (i in seq_along(terms_cycle_through)) {
  
  # exclude feature i
  currentFeatures <- terms_cycle_through[-i]
  
  # assemble regression formula
  regressionFormula <- as.formula(paste("expression ~ ", paste(currentFeatures, collapse="+")))
  
  # fit model
  currentModel <- lm(formula = regressionFormula, data = all_model_data)
  
  # Calculate difference in R2 value with and without the variable
  term_result <- data.frame(term = terms_cycle_through[i], r2 = (full_model_r2 - broom::glance(currentModel)$adj.r.squared))
  outsideLoopTermResult <- rbind(outsideLoopTermResult, term_result)
}

# View results
View(outsideLoopTermResult)


# Create function to test individual TE terms ----------------------------------

# Create a function that runs the same model 
modelDistance <- function(variableName, umrVariableName, simplifiedName, df) {
  # Create model statement
  modelStatement <- as.formula(paste0("expression ~ ", variableName, " + ", umrVariableName, " + mean_expression +  mean_phylop + mean_dnds + age + tissue"))
  print(modelStatement)
  
  # Run model
  runModel <- lm(modelStatement, data = df)
  
  # Collect model metrics
  individualFullModelStats <- data.frame(adj_r2 = broom::glance(runModel)$adj.r.squared,
                                         pvalue = broom::glance(runModel)$p.value,
                                         model = paste0(simplifiedName, " Model"))
  
  # Collect term effects
  fullModelnoUMREffects <- broom::tidy(runModel)
  fullModelnoUMREffects$model <- rep(paste0(simplifiedName, " Model"), rep(nrow(fullModelnoUMREffects)))
  
  # Save both as list
  return(list(individualFullModelStats, fullModelnoUMREffects))
}


# Use the function to run the data ---------------------------------------------

m1 <- modelDistance("scaledTEBP_1_500", "umrCount_1_500", "1-500", all_model_data)
m2 <- modelDistance("scaledTEBP_501_1000", "umrCount_501_1000", "501_1000", all_model_data)
m3 <- modelDistance("scaledTEBP_1001_1500", "umrCount_1001_1500", "1001_1500", all_model_data)
m4 <- modelDistance("scaledTEBP_1501_3000", "umrCount_1501_3000", "1501_3000", all_model_data)
m5 <- modelDistance("scaledTEBP_3001_4500", "umrCount_3001_4500", "3001_4500", all_model_data)
m6 <- modelDistance("scaledTEBP_4501_6000", "umrCount_4501_6000", "4501_6000", all_model_data)
m7 <- modelDistance("scaledTEBP_6001_7500", "umrCount_6001_7500", "6001_7500", all_model_data)
m8 <- modelDistance("scaledTEBP_7501_9000", "umrCount_7501_9000", "7501_9000", all_model_data)
m9 <- modelDistance("scaledTEBP_9001_10500", "umrCount_9001_10500", "9001_10500", all_model_data)
m10 <- modelDistance("scaledTEBP_10501_12000", "umrCount_10501_12000", "10501_12000", all_model_data)
m11 <- modelDistance("scaledTEBP_12001_13500", "umrCount_12001_13500", "12001_13500", all_model_data)
m12 <- modelDistance("scaledTEBP_13501_15000", "umrCount_13501_15000", "13501_15000", all_model_data)


# Combine all model results ----------------------------------------------------

# Partial effect estimates
allModelSummaries <- rbind(jointModelnoUMREffects, 
                           m1[[2]], m2[[2]], m3[[2]], m4[[2]], m5[[2]], m6[[2]], 
                           m7[[2]], m8[[2]], m9[[2]], m10[[2]], m11[[2]], m12[[2]])

# Full model stats
allFullModelStats <- rbind(jointModelStats, 
                           m1[[1]], m2[[1]], m3[[1]], m4[[1]], m5[[1]], m6[[1]], 
                           m7[[1]], m8[[1]], m9[[1]], m10[[1]], m11[[1]], m12[[1]])

# Export to file
write.csv(allModelSummaries, 
          "/workdir/mbb262/model_data/inbred_B73Only_model_withUMRs_effectStatistics_single_joint.csv",
          row.names = FALSE,
          quote = FALSE)
write.csv(allFullModelStats, 
          "/workdir/mbb262/model_data/inbred_B73Only_model_withUMRs_modelFit_single_joint.csv",
          row.names = FALSE,
          quote = FALSE)


# Plot model results -----------------------------------------------------------

# Subet file to just the distance terms
allDistanceTerms <- c("scaledTEBP_1_500", "scaledTEBP_501_1000", "scaledTEBP_1001_1500", 
                      "scaledTEBP_1501_3000", "scaledTEBP_3001_4500", "scaledTEBP_4501_6000", 
                      "scaledTEBP_6001_7500", "scaledTEBP_7501_9000", "scaledTEBP_9001_10500", 
                      "scaledTEBP_10501_12000", "scaledTEBP_12001_13500", "scaledTEBP_13501_15000",
                      "umrCount_1_500", "umrCount_501_1000", "umrCount_1001_1500", 
                      "umrCount_1501_3000", "umrCount_3001_4500", "umrCount_4501_6000", 
                      "umrCount_6001_7500", "umrCount_7501_9000", "umrCount_9001_10500", 
                      "umrCount_10501_12000", "umrCount_12001_13500", "umrCount_13501_15000")
sub_allModelSummaries <- allModelSummaries %>% 
  filter(term %in% allDistanceTerms)

# Create some labels for plotting
sub_allModelSummaries$modelType <- gsub("^[0-9].*", "Single Term Model", sub_allModelSummaries$model)
sub_allModelSummaries$Model <- paste0(gsub("_.*", "", sub_allModelSummaries$term), " ", sub_allModelSummaries$modelType)
sub_allModelSummaries$Model <- gsub("scaled", "", sub_allModelSummaries$Model)
sub_allModelSummaries$Model <- gsub("umrCount", "UMR", sub_allModelSummaries$Model)

sub_allModelSummaries$Significant <- sub_allModelSummaries$p.value <= 0.05
sub_allModelSummaries$Term <- gsub("scaledTEBP_", "", sub_allModelSummaries$term)
sub_allModelSummaries$Term <- gsub("umrCount_", "", sub_allModelSummaries$Term)
sub_allModelSummaries$Term <- gsub("_", "-", sub_allModelSummaries$Term)

# Change ordering of factors
sub_allModelSummaries$Term <- as.character(sub_allModelSummaries$Term)
sub_allModelSummaries$Term <- factor(sub_allModelSummaries$Term, 
                                     levels=c("1-500","501-1000", "1001-1500", 
                                              "1501-3000", "3001-4500", "4501-6000",
                                              "6001-7500", "7501-9000", "9001-10500", 
                                              "10501-12000", "12001-13500", "13501-15000"))


# Make plots -------------------------------------------------------------------

# Working example: "R**2*'='*0.9*', '*p<2.2e-16"
plotText <- paste0("'Joint Model '*R**2*'='*", 
                   allFullModelStats %>% filter(model == "Joint Term Model") %>% select(adj_r2) %>% round(digits = 3),
                   "*', '*p<2.2e-16")
r2_text_x = 1
plotText

# Plot plot
plot1 <- ggplot(sub_allModelSummaries, aes(y = estimate, x = Term, group = Model)) +
  geom_line(aes(color=Model)) +
  geom_point(aes(shape=Significant)) +
  xlab("Distance from Gene") +
  ylab("Effect") +
  theme_classic() +
  annotate(geom="text", label = plotText, x = 1.5, y = 0.047, size = annotate_size, parse = TRUE, hjust = 0) +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.box="vertical",
        legend.text = element_text(size = axis_text_size),
        legend.title=element_text(size = legend_text_size),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette)
plot1


# Save plot --------------------------------------------------------------------

ggsave("/home/mbb262/git_projects/te_ase_nam/figs/inbred_B73Only_acrossFamilies_model_effect_joint_single_term_withUMRs.png",
       plot = plot1,
       width = 7.5,
       height = 5.5, 
       units = "in",
       dpi = "retina")


