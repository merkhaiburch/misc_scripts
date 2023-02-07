# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-06-06
# Updated... 2022-06-06
#
# Description: Model 1 (RF with no GO, only nuisance/mapping terms)
# Using pleiotropy scores from number of unique traits
# combined with ONLY nuisance/mapping terms in a random forest model
# 
# Running 8 models, 4: Observed data, 4: permuted data
# Each of the 4 models comes from the 4 population-trait categories
#
# Model (the same one used to calculate adjusted scores):
# pleiotropy score ~ average LD + interval size + input SNP count
# ---------------------------------------------------------------

# tutorials followed for ranger::ranger
# https://brunaw.com/slides/rladies-dublin/RF/intro-to-rf.html#31
# https://uc-r.github.io/random_forests

# Load in packages
library(dplyr)
library(tidyr)
library(ranger)
library(ggplot2)

# set seed
set.seed(2022)

# read in data aggregatied by aggregate_interval.data.R
data_dir <- "~/git_projects/pleiotropy/data/"
all_interval_data <- data.table::fread(paste0(data_dir, "/all_interval_data_unique_filter.csv"))

# Replace NAs with -1 for GERP columns
all_interval_data$np_noperm_mean_gerp[is.na(all_interval_data$np_noperm_mean_gerp)] <- -1
all_interval_data$gp_noperm_mean_gerp[is.na(all_interval_data$gp_noperm_mean_gerp)] <- -1
all_interval_data$gm_noperm_mean_gerp[is.na(all_interval_data$gm_noperm_mean_gerp)] <- -1
all_interval_data$ge_noperm_mean_gerp[is.na(all_interval_data$ge_noperm_mean_gerp)] <- -1


# ---------------------------------
#           Format data
# ---------------------------------

# Melt data for permuted columns
all_interval_data_melt_permuted <- all_interval_data %>% 
  select("rr_id", 
         "range_type", 
         "seqid",
         "interval_size",
         "average_r2_ld_nam",
         "average_r2_ld_282",
         "snp_count_nam_all",
         "snp_count_goodman_all",
         "npp_unique_perm_1","npp_unique_perm_2","npp_unique_perm_3","npp_unique_perm_4","npp_unique_perm_5",
         "npp_unique_perm_6","npp_unique_perm_7","npp_unique_perm_8","npp_unique_perm_9","npp_unique_perm_10",
         "gpp_unique_perm_1","gpp_unique_perm_2","gpp_unique_perm_3","gpp_unique_perm_4","gpp_unique_perm_5",
         "gpp_unique_perm_6","gpp_unique_perm_7","gpp_unique_perm_8","gpp_unique_perm_9","gpp_unique_perm_10",
         "gpm_perm_1","gpm_perm_2","gpm_perm_3","gpm_perm_4","gpm_perm_5",
         "gpm_perm_6","gpm_perm_7","gpm_perm_8","gpm_perm_9","gpm_perm_10",
         "gpe_perm_1","gpe_perm_2","gpe_perm_3","gpe_perm_4","gpe_perm_5")


# Melt individual permuted datasets, then merge 
nam_field <- tidyr::pivot_longer(all_interval_data_melt_permuted,
                                 cols = npp_unique_perm_1:npp_unique_perm_10,
                                 names_to = "pop_trait_perm_no_perm_level",
                                 values_to = "pleiotropy_npp")
goodman_field <- tidyr::pivot_longer(all_interval_data_melt_permuted,
                                     cols = gpp_unique_perm_1:gpp_unique_perm_10,
                                     names_to = "pop_trait_perm_no_perm_level",
                                     values_to = "pleiotropy_gpp")
goodman_mass <- tidyr::pivot_longer(all_interval_data_melt_permuted,
                                    cols = gpm_perm_1:gpm_perm_10,
                                    names_to = "pop_trait_perm_no_perm_level",
                                    values_to = "pleiotropy_gpm")

# Add a column with information about which permutation
nam_field$perm_number <- gsub(".*_", "", nam_field$pop_trait_perm_no_perm_level)
goodman_field$perm_number <- gsub(".*_", "", goodman_field$pop_trait_perm_no_perm_level)
goodman_mass$perm_number <- gsub(".*_", "", goodman_mass$pop_trait_perm_no_perm_level)

# Merge
permuted_merged <- Reduce(function(x, y) merge(x, y, by = c("rr_id", "perm_number", "range_type")), 
                          list(nam_field[,c(1:8,35:36)], 
                               goodman_field[,c(1:2,35:36)], 
                               goodman_mass[,c(1:2,35:36)]))
dim(permuted_merged) # should be: 754900     12
colnames(permuted_merged)

# Keep this one separate because it only has 5 permutations
goodman_expression <- tidyr::pivot_longer(all_interval_data_melt_permuted,
                                          cols = gpe_perm_1:gpe_perm_5,
                                          names_to = "pop_trait_perm_no_perm_level",
                                          values_to = "pleiotropy_gpe")



# --------------------------------------------
# Partition data into Observed Data and Permuted sets
# and separate population-trait categories
# --------------------------------------------

# ------------------------------
# nam Observed Data
# ------------------------------

np_obs <- all_interval_data %>% select("rr_id", 
                                       "seqid", 
                                       "average_r2_ld_nam", 
                                       "snp_count_nam_all", 
                                       "interval_size",
                                       "nam_filtered_all_true_uniqueCount")
colnames(np_obs)
dim(np_obs) # 75490    6

# ------------------------------
# NAM Permuted
# ------------------------------

np_permuted <- permuted_merged %>% select("rr_id", 
                                          "seqid", 
                                          "average_r2_ld_nam", 
                                          "snp_count_nam_all", 
                                          "interval_size",
                                          "pleiotropy_npp")
colnames(np_permuted)
dim(np_permuted) # 754,900     6


# ------------------------------
# goodman Field Observed Data
# ------------------------------

gp_obs <- all_interval_data %>% select("rr_id", 
                                       "seqid", 
                                       "average_r2_ld_282", 
                                       "snp_count_goodman_all", 
                                       "interval_size",
                                       "goodman_filtered_physiological_true_uniqueCount")
colnames(gp_obs)
dim(gp_obs) # 75490    6


# ------------------------------
# Goodman Field Permuted
# ------------------------------

gp_permuted <- permuted_merged %>% select("rr_id", 
                                          "seqid", 
                                          "average_r2_ld_282", 
                                          "snp_count_goodman_all", 
                                          "interval_size",
                                          "pleiotropy_gpp")
colnames(gp_permuted)
dim(gp_permuted) # 754,900     6


# ------------------------------
# goodman Mass Features Observed Data
# ------------------------------

gm_obs <- all_interval_data %>% select("rr_id", 
                                       "seqid", 
                                       "average_r2_ld_282", 
                                       "snp_count_goodman_all", 
                                       "interval_size",
                                       "goodman_filtered_metabolite_true_uniqueCount")
colnames(gm_obs)
dim(gm_obs)# 75490    6


# ------------------------------
# Goodman Mass Features Permuted
# ------------------------------

gm_permuted <- permuted_merged %>% select("rr_id", 
                                          "seqid", 
                                          "average_r2_ld_282", 
                                          "snp_count_goodman_all", 
                                          "interval_size",
                                          "pleiotropy_gpm")
colnames(gm_permuted)
dim(gm_permuted) # 754,900     6


# ------------------------------
# goodman expression Observed Data
# ------------------------------

ge_obs <- all_interval_data %>% select("rr_id", 
                                       "seqid", 
                                       "average_r2_ld_282", 
                                       "snp_count_goodman_all", 
                                       "interval_size",
                                       "goodman_filtered_expression_uniqueCount")
colnames(ge_obs)
dim(ge_obs) # 75490    6


# ------------------------------
# Goodman Expression Permuted
# ------------------------------

ge_permuted <- goodman_expression %>% select("rr_id", 
                                             "seqid", 
                                             "average_r2_ld_282", 
                                             "snp_count_goodman_all", 
                                             "interval_size",
                                             "pleiotropy_gpe")
colnames(ge_permuted)
dim(ge_permuted) # 377,450     6


# ------------------------------------------------------------------------------
# Split training and testing data
# FILTER BASED ON CHROMOSOME, 1:9 TRAIN, 10 TEST
# Doing this through the case weights and holdout method in ranger::ranger()
# ------------------------------------------------------------------------------

# nam Field - Observed Data
np_train_obs <- np_obs %>% arrange(seqid)
np_test_obs <- np_obs %>% filter(seqid == 10)
np_case_obs <- c(rep(1, nrow(np_obs %>% filter(seqid !=10))),
                 rep(0, nrow(np_obs %>% filter(seqid == 10))))

# nam Field - Permuted
np_train_permuted <- np_permuted %>% arrange(seqid)
np_test_permuted <- np_permuted %>% filter(seqid == 10)
np_case_permuted <- c(rep(1, nrow(np_permuted %>% filter(seqid !=10))),
                      rep(0, nrow(np_permuted %>% filter(seqid == 10))))

# 282 Field - Observed Data
gp_train_obs <- gp_obs %>% arrange(seqid)
gp_test_obs <- gp_obs %>% filter(seqid == 10)
gp_case_obs <- c(rep(1, nrow(gp_obs %>% filter(seqid !=10))),
                 rep(0, nrow(gp_obs %>% filter(seqid == 10))))

# 282 Field - Permuted
gp_train_permuted <- gp_permuted %>% arrange(seqid)
gp_test_permuted <- gp_permuted %>% filter(seqid == 10)
gp_case_permuted <- c(rep(1, nrow(gp_permuted %>% filter(seqid !=10))),
                      rep(0, nrow(gp_permuted %>% filter(seqid == 10))))

# 282 Mass Features - Observed Data
gm_test_obs <- gm_obs %>% arrange(seqid)
gm_test_obs <- gm_obs %>% filter(seqid == 10)
gm_case_obs <- c(rep(1, nrow(gm_obs %>% filter(seqid !=10))),
                 rep(0, nrow(gm_obs %>% filter(seqid == 10))))

# 282 Mass Features - Permuted
gm_train_permuted <- gm_permuted %>% arrange(seqid)
gm_test_permuted <- gm_permuted %>% filter(seqid == 10)
gm_case_permuted <- c(rep(1, nrow(gm_permuted %>% filter(seqid !=10))),
                      rep(0, nrow(gm_permuted %>% filter(seqid == 10))))

# 282 expression - Observed Data
ge_train_obs <- ge_obs %>% arrange(seqid)
ge_test_obs <- ge_obs %>% filter(seqid == 10)
ge_case_obs <- c(rep(1, nrow(ge_obs %>% filter(seqid !=10))),
                 rep(0, nrow(ge_obs %>% filter(seqid == 10))))

# 282 expression - Permuted
ge_train_permuted <- ge_permuted %>% arrange(seqid)
ge_test_permuted <- ge_permuted %>% filter(seqid == 10)
ge_case_permuted <- c(rep(1, nrow(ge_permuted %>% filter(seqid !=10))),
                      rep(0, nrow(ge_permuted %>% filter(seqid == 10))))


# ------------------------------------------------
# Format names to look nicer for RF model outputs
# ------------------------------------------------

# load in table with formatted names
pretty_names <- data.table::fread("~/git_projects/pleiotropy/data/pretty_names_for_rf.csv")


# ---------------------------------------------------------------------------------
#                       Run Random Forest Models
# ---------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
#                               nam Field 
# --------------------------------------------------------------------------------

# Observed Data model ----------------------------------
nam_field_obs <- ranger::ranger(nam_filtered_all_true_uniqueCount ~ ., 
                                importance = 'impurity',
                                case.weights = np_case_obs,
                                holdout = TRUE,
                                data = np_obs[,-c(1,2)],
                                num.trees = 500)
nam_field_obs

# prep for importance plot
nam_field_obs_importance <- as.data.frame(nam_field_obs$variable.importance)
nam_field_obs_importance$type <- rownames(nam_field_obs_importance)
colnames(nam_field_obs_importance) <- c("Importance", "Variable_Type")
nam_field_obs_importance$Importance <- nam_field_obs_importance$Importance/max(nam_field_obs_importance$Importance)
nam_field_obs_importance <- merge(x = nam_field_obs_importance, y = pretty_names, 
                                  by.x = "Variable_Type", by.y = "old", all.x = TRUE)

# Permuted model-------------------------------------
nam_field_perm <- ranger::ranger(pleiotropy_npp ~ ., 
                                 importance = 'impurity',
                                 case.weights = np_case_permuted,
                                 holdout = TRUE,
                                 data = np_permuted[,-c(1,2)],
                                 num.trees = 500)
nam_field_perm

# prep for importance plot
nam_field_perm_importance <- as.data.frame(nam_field_perm$variable.importance)
nam_field_perm_importance$type <- rownames(nam_field_perm_importance)
colnames(nam_field_perm_importance) <- c("Importance", "Variable_Type")
nam_field_perm_importance$Importance <- nam_field_perm_importance$Importance/max(nam_field_perm_importance$Importance)
nam_field_perm_importance <- merge(x = nam_field_perm_importance, y = pretty_names, 
                                   by.x = "Variable_Type", by.y = "old", all.x = TRUE)


# --------------------------------------------------------------------------------
#                               282 Field
# --------------------------------------------------------------------------------

# 282 Field Observed Data-----------------------------
goodman_field_obs <- ranger::ranger(goodman_filtered_physiological_true_uniqueCount ~ ., 
                                    importance = 'impurity',
                                    case.weights = gp_case_obs,
                                    holdout = TRUE,
                                    data = gp_obs[,-c(1,2)],
                                    num.trees = 500)
goodman_field_obs

# prep for importance plot
goodman_field_obs_importance <- as.data.frame(goodman_field_obs$variable.importance)
goodman_field_obs_importance$type <- rownames(goodman_field_obs_importance)
colnames(goodman_field_obs_importance) <- c("Importance", "Variable_Type")
goodman_field_obs_importance$Importance <- goodman_field_obs_importance$Importance/max(goodman_field_obs_importance$Importance)
goodman_field_obs_importance <- merge(x = goodman_field_obs_importance, y = pretty_names, 
                                      by.x = "Variable_Type", by.y = "old", all.x = TRUE)


# 282 Field Permuted----------------------------------
goodman_field_perm <- ranger::ranger(pleiotropy_gpp ~ ., 
                                     importance = 'impurity',
                                     case.weights = gp_case_permuted,
                                     holdout = TRUE,
                                     data = gp_permuted[,-c(1,2)],
                                     num.trees = 500)
goodman_field_perm

# prep for importance plot
goodman_field_perm_importance <- as.data.frame(goodman_field_perm$variable.importance)
goodman_field_perm_importance$type <- rownames(goodman_field_perm_importance)
colnames(goodman_field_perm_importance) <- c("Importance", "Variable_Type")
goodman_field_perm_importance$Importance <- goodman_field_perm_importance$Importance/max(goodman_field_perm_importance$Importance)
goodman_field_perm_importance <- merge(x = goodman_field_perm_importance, y = pretty_names, 
                                       by.x = "Variable_Type", by.y = "old", all.x = TRUE)


# --------------------------------------------------------------------------------
#                               282 Mass Features
# --------------------------------------------------------------------------------

# Observed Data model-----------------------------------
goodman_mass_obs <- ranger::ranger(goodman_filtered_metabolite_true_uniqueCount ~ ., 
                                   importance = 'impurity',
                                   case.weights = gm_case_obs,
                                   holdout = TRUE,
                                   data = gm_obs[,-c(1,2)],
                                   num.trees = 500)
goodman_mass_obs

# prep for importance plot
goodman_mass_obs_importance <- as.data.frame(goodman_mass_obs$variable.importance)
goodman_mass_obs_importance$type <- rownames(goodman_mass_obs_importance)
colnames(goodman_mass_obs_importance) <- c("Importance", "Variable_Type")
goodman_mass_obs_importance$Importance <- goodman_mass_obs_importance$Importance/max(goodman_mass_obs_importance$Importance)
goodman_mass_obs_importance <- merge(x = goodman_mass_obs_importance, y = pretty_names, 
                                     by.x = "Variable_Type", by.y = "old", all.x = TRUE)

# Permuted model-----------------------------------------
goodman_mass_perm <- ranger::ranger(pleiotropy_gpm ~ ., 
                                    importance = 'impurity',
                                    case.weights = gm_case_permuted,
                                    holdout = TRUE,
                                    data = gm_permuted[,-c(1,2)],
                                    num.trees = 500)
goodman_mass_perm

# prep for importance plot
goodman_mass_perm_importance <- as.data.frame(goodman_mass_perm$variable.importance)
goodman_mass_perm_importance$type <- rownames(goodman_mass_perm_importance)
colnames(goodman_mass_perm_importance) <- c("Importance", "Variable_Type")
goodman_mass_perm_importance$Importance <- goodman_mass_perm_importance$Importance/max(goodman_mass_perm_importance$Importance)
goodman_mass_perm_importance <- merge(x = goodman_mass_perm_importance, y = pretty_names, 
                                      by.x = "Variable_Type", by.y = "old", all.x = TRUE)


# --------------------------------------------------------------------------------
#                                 282 expression
# --------------------------------------------------------------------------------

# Observed Data model----------------------------------------------
goodman_expr_obs <- ranger::ranger(goodman_filtered_expression_uniqueCount ~ ., 
                                   importance = 'impurity',
                                   case.weights = ge_case_obs,
                                   holdout = TRUE,
                                   data = ge_obs[,-c(1,2)],
                                   num.trees = 500)
goodman_expr_obs

# prep for importance plot
goodman_expr_obs_importance <- as.data.frame(goodman_expr_obs$variable.importance)
goodman_expr_obs_importance$type <- rownames(goodman_expr_obs_importance)
colnames(goodman_expr_obs_importance) <- c("Importance", "Variable_Type")
goodman_expr_obs_importance$Importance <- goodman_expr_obs_importance$Importance/max(goodman_expr_obs_importance$Importance)
goodman_expr_obs_importance <- merge(x = goodman_expr_obs_importance, y = pretty_names, 
                                     by.x = "Variable_Type", by.y = "old", all.x = TRUE)

# Permuted model ------------------------------------------------------
goodman_expr_perm <- ranger::ranger(pleiotropy_gpe ~ ., 
                                    importance = 'impurity',
                                    case.weights = ge_case_permuted,
                                    holdout = TRUE,
                                    data = ge_permuted[,-c(1,2)],
                                    num.trees = 500)
goodman_expr_perm

# prep for importance plot
goodman_expr_perm_importance <- as.data.frame(goodman_expr_perm$variable.importance)
goodman_expr_perm_importance$type <- rownames(goodman_expr_perm_importance)
colnames(goodman_expr_perm_importance) <- c("Importance", "Variable_Type")
goodman_expr_perm_importance$Importance <- goodman_expr_perm_importance$Importance/max(goodman_expr_perm_importance$Importance)
goodman_expr_perm_importance <- merge(x = goodman_expr_perm_importance, y = pretty_names, 
                                      by.x = "Variable_Type", by.y = "old", all.x = TRUE)


# --------------------------------------------------------------------------------
# Make a table summarizing results
# --------------------------------------------------------------------------------

nam_field_obs # nam field Observed data
nam_field_perm # nam physiological/field permuted
goodman_field_obs
goodman_field_perm
goodman_mass_obs
goodman_mass_perm
goodman_expr_obs
goodman_expr_perm

# Assemble model results
temp <- rbind(nam_field_obs[c("prediction.error", "r.squared")] %>% data.frame(), 
              nam_field_perm[c("prediction.error", "r.squared")] %>% data.frame(),
              goodman_field_obs[c("prediction.error", "r.squared")] %>% data.frame(),
              goodman_field_perm[c("prediction.error", "r.squared")] %>% data.frame(),
              goodman_mass_obs[c("prediction.error", "r.squared")] %>% data.frame(),
              goodman_mass_perm[c("prediction.error", "r.squared")] %>% data.frame(),
              goodman_expr_obs[c("prediction.error", "r.squared")] %>% data.frame(),
              goodman_expr_perm[c("prediction.error", "r.squared")] %>% data.frame())
temp$model_id <- c("NAM Field Observed Data", 
                   "NAM Field Permuted Data",
                   "GAP Field Observed Data",
                   "GAP Field Permuted Data",
                   "GAP Mass Features Observed Data",
                   "GAP Mass Features Permuted Data",
                   "GAP Expression Observed Data",
                   "GAP Expression Permuted Data")
# rearrange and make nicer names
temp <- temp[,c(3,1,2)]
colnames(temp) <- c("model", "prediction_error", "r_squared")

# Export
data.table::fwrite(temp, "~/git_projects/pleiotropy/data/model1_RF_prediction_accuracy.csv")


# -------------------------
# Set plot size parameters
# -------------------------

axis_text_size <- 15
title_text_size <- 16
tag_size <- 16
legend_text_size <- 14
legend_shape_size <- 0.95
r2_text_x <- 0
annotate_size <- 3.5

# The color blind palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# --------------------------------------------------------------------------------
#                         Feature importance plots 
# --------------------------------------------------------------------------------

# NAM Phys Observed Data -----------------------------------------------
a <- nam_field_obs_importance %>% 
  dplyr::arrange(desc(Importance)) %>%
  ggplot(aes(x= reorder(Variable, Importance), y = Importance, fill = Type)) +
  geom_col() +
  coord_flip() +
  ggtitle("NAM Field Observed") +
  xlab("Variable Type") +
  ylab("Relative Importance") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9","#E69F00")) +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
  labs(tag = "A")


# nam phys Permuted -------------------------------------------------
b <- nam_field_perm_importance %>% 
  dplyr::arrange(desc(Importance)) %>%
  ggplot(aes(x= reorder(Variable, Importance), y = Importance, fill = Type)) +
  geom_col() +
  coord_flip() +
  ggtitle("NAM Field Permuted") +
  xlab("Variable Type") +
  ylab("Relative Importance") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9","#E69F00")) +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
  labs(tag = "A")

# goodman Field Observed Data ----------------------------------
c <- goodman_field_obs_importance %>% 
  dplyr::arrange(desc(Importance)) %>%
  ggplot(aes(x= reorder(Variable, Importance), y = Importance, fill = Type)) +
  geom_col() +
  coord_flip() +
  ggtitle("GAP Field Observed") +
  xlab("Variable Type") +
  ylab("Relative Importance")+
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9","#E69F00")) +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
  labs(tag = "B")

# goodman Field Permuted ------------------------------------
d <- goodman_field_perm_importance %>% 
  dplyr::arrange(desc(Importance)) %>%
  ggplot(aes(x= reorder(Variable, Importance), y = Importance, fill = Type)) +
  geom_col() +
  coord_flip() +
  ggtitle("GAP Field Permuted") +
  xlab("Variable Type") +
  ylab("Relative Importance")+
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9","#E69F00")) +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
  labs(tag = "B")

# goodman Mass Features Observed Data -------------------------------------
e <- goodman_mass_obs_importance %>% 
  dplyr::arrange(desc(Importance)) %>%
  ggplot(aes(x= reorder(Variable, Importance), y = Importance, fill = Type)) +
  geom_col() +
  coord_flip() +
  ggtitle("GAP Mass Features Observed") +
  xlab("Variable Type") +
  ylab("Relative Importance") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9","#E69F00")) +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
  labs(tag = "C")

# goodman Mass Features Permuted --------------------------------------
f <- goodman_mass_perm_importance %>% 
  dplyr::arrange(desc(Importance)) %>%
  ggplot(aes(x= reorder(Variable, Importance), y = Importance, fill = Type)) +
  geom_col() +
  coord_flip() +
  ggtitle("GAP Mass Features Permuted") +
  xlab("Variable Type") +
  ylab("Relative Importance") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9","#E69F00")) +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
  labs(tag = "C")

# goodman expression Observed Data ------------------------------------
g <- goodman_expr_obs_importance %>% 
  dplyr::arrange(desc(Importance)) %>%
  ggplot(aes(x= reorder(Variable, Importance), y = Importance, fill = Type)) +
  geom_col() +
  coord_flip() +
  ggtitle("GAP Expression Observed") +
  xlab("Variable Type") +
  ylab("Relative Importance") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9")) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  labs(tag = "D")

# goodman expression Permuted --------------------------------------
h <- goodman_expr_perm_importance %>% 
  dplyr::arrange(desc(Importance)) %>%
  ggplot(aes(x= reorder(Variable, Importance), y = Importance, fill = Type)) +
  geom_col() +
  coord_flip() +
  ggtitle("GAP Expression Permuted") +
  xlab("Variable Type") +
  ylab("Relative Importance") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9")) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  labs(tag = "D")


# test combined plot ----------------------------------------------
# ggpubr::ggarrange(a,b,c,d,
#                   e,f,g,h,
#                   nrow = 4, ncol = 2,
#                   common.legend = TRUE, legend = "bottom", align = "v")

# Export to file---------------------------------------------------
ggsave("~/git_projects/pleiotropy/images/model1_random_forest_importance_observed_vs_permuted.png",
       plot = ggpubr::ggarrange(a,b,c,d,e,f,g,h,
                                nrow = 4, ncol = 2, common.legend = TRUE, legend = "none", align = "v"),
       width = 15,
       height = 17, 
       units = "in",
       dpi = "retina")

# Observed data
ggsave("~/git_projects/pleiotropy/images/model1_random_forest_importance_observed.png",
       plot = ggpubr::ggarrange(a,c,e,g, nrow = 2, ncol = 2, 
                                common.legend = TRUE, legend = "none", align = "v"),
       width = 17,
       height = 13, 
       units = "in",
       dpi = "retina")

# Permuted data
ggsave("~/git_projects/pleiotropy/images/model1_random_forest_importance_permuted.png",
       plot = ggpubr::ggarrange(b,d,f,h,nrow = 2, ncol = 2, 
                                common.legend = TRUE, legend = "none", align = "v"),
       width = 17,
       height = 13, 
       units = "in",
       dpi = "retina")


# --------------------------------------------------------------------------------
# Plot & combine all model performance plots 
# --------------------------------------------------------------------------------

# nam phys Observed Data --------------------------------------------------------------
np_pred_obs <- predict(object = nam_field_obs, data = np_test_obs)
np_test_obs$np_pred_obs <- np_pred_obs$predictions

test2 <- paste("~R^2==~", round(cor(np_test_obs$nam_filtered_all_true_uniqueCount, np_test_obs$np_pred_obs)^2, digits = 3))

a_pred <- ggplot(np_test_obs, aes(x = nam_filtered_all_true_uniqueCount, y = np_pred_obs)) +
  geom_point(colour = "#999999") +
  geom_smooth(method = lm, color = "black", linetype = "solid") + 
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  labs(x = "Observed Pleiotropy", 
       y = "Predicted Pleiotropy", 
       title = "NAM Field Observed") +
  annotate(geom="text", label = test2, x = r2_text_x, y = 60, size = annotate_size, parse = TRUE, hjust = 0) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  ylim(0,65) +
  xlim(0,65) +
  labs(tag = "A")


# nam phys Permuted-----------------------------------------------------------------
np_pred_permuted <- predict(object = nam_field_perm, data = np_test_permuted)
np_test_permuted$np_pred_permuted <- np_pred_permuted$predictions

test2 <- paste("~R^2==~", 
               round(cor(np_test_permuted$pleiotropy_npp, np_test_permuted$np_pred_permuted)^2, digits = 3))

b_pred <- ggplot(np_test_permuted, aes(x = pleiotropy_npp, y = np_pred_permuted)) +
  geom_point(colour = "#999999") +
  geom_smooth(method = lm, color = "black") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Permuted Pleiotropy", 
       y = "Predicted Pleiotropy", 
       title = "NAM Field Permuted") +
  annotate(geom="text", label = test2, x = r2_text_x, y = 60, size = annotate_size, parse = TRUE, hjust = 0) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  ylim(0,65) +
  xlim(0,65) +
  labs(tag = "A")

# goodman Field Observed Data--------------------------------------------------
gp_pred_obs <- predict(object = goodman_field_obs, data = gp_test_obs)
gp_test_obs$gp_pred_obs <- gp_pred_obs$predictions

test2 <- paste("~R^2==~", 
               round(cor(gp_test_obs$goodman_filtered_physiological_true_uniqueCount, gp_test_obs$gp_pred_obs)^2, digits = 3))

c_pred <- ggplot(gp_test_obs, aes(x = goodman_filtered_physiological_true_uniqueCount, y = gp_pred_obs)) +
  geom_point(colour = "#999999") +
  geom_smooth(method = lm, color = "black") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Observed Pleiotropy", 
       y = "Predicted Pleiotropy", 
       title = "GAP Field Observed") +
  annotate(geom="text", label = test2, x = r2_text_x, y = 60, size = annotate_size, parse = TRUE, hjust = 0) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  ylim(0,65) +
  xlim(0,65) +
  labs(tag = "B")

# goodman Field Permuted---------------------------------------------------
gp_pred_permuted <- predict(object = goodman_field_perm, data = gp_test_permuted)
gp_test_permuted$gp_pred_permuted <- gp_pred_permuted$predictions

test2 <- paste("~R^2==~", 
               round(cor(gp_test_permuted$pleiotropy_gpp, gp_test_permuted$gp_pred_permuted)^2, digits = 3))

d_pred <- ggplot(gp_test_permuted, aes(x = pleiotropy_gpp, y = gp_pred_permuted)) +
  geom_point(colour = "#999999") +
  geom_smooth(method = lm, color = "black") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Permuted Pleiotropy", 
       y = "Predicted Pleiotropy", 
       title = "GAP Field Permuted") +
  annotate(geom="text", label = test2, x = r2_text_x, y = 60, size = annotate_size, parse = TRUE, hjust = 0) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  ylim(0,65) +
  xlim(0,65) +
  labs(tag = "B")

# goodman Mass Features Observed Data----------------------------------------------------
# Predictions
gm_pred_obs <- predict(object = goodman_mass_obs, data = gm_test_obs)
gm_test_obs$gm_pred_obs <- gm_pred_obs$predictions

test2 <- paste("~R^2==~", 
               round(cor(gm_test_obs$goodman_filtered_metabolite_true_uniqueCount, gm_test_obs$gm_pred_obs)^2, digits = 3))

e_pred <- ggplot(gm_test_obs, aes(x = goodman_filtered_metabolite_true_uniqueCount, y = gm_pred_obs)) +
  geom_point(colour = "#999999") +
  geom_smooth(method = lm, color = "black") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Observed Pleiotropy", 
       y = "Predicted Pleiotropy", 
       title = "GAP Mass Features Observed") +
  annotate(geom="text", label = test2, x = r2_text_x, y = 1500, size = annotate_size, parse = TRUE, hjust = 0) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  ylim(0,1850) +
  xlim(0,1850) +
  labs(tag = "C")

# goodman Mass Features Permuted-----------------------------------------------------
gm_pred_permuted <- predict(object = goodman_mass_perm, data = gm_test_permuted)
gm_test_permuted$gm_pred_permuted <- gm_pred_permuted$predictions

test2 <- paste("~R^2==~", 
               round(cor(gm_test_permuted$pleiotropy_gpm, gm_test_permuted$gm_pred_permuted)^2, digits = 3))

f_pred <- ggplot(gm_test_permuted, aes(x = pleiotropy_gpm, y = gm_pred_permuted)) +
  geom_point(colour = "#999999") +
  geom_smooth(method = lm, color = "black") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Permuted Pleiotropy", 
       y = "Predicted Pleiotropy", 
       title = "GAP Mass Features Permuted") +
  annotate(geom="text", label = test2, x = r2_text_x, y = 1500, size = annotate_size, parse = TRUE, hjust = 0) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  ylim(0,1850) +
  xlim(0,1850) +
  labs(tag = "C")

# goodman expression Observed Data---------------------------------------------------
ge_pred_obs <- predict(object = goodman_expr_obs, data = ge_test_obs)
ge_test_obs$ge_pred_obs <- ge_pred_obs$predictions

test2 <- paste("~R^2==~", 
               round(cor(ge_test_obs$goodman_filtered_expression_uniqueCount, ge_test_obs$ge_pred_obs)^2, digits = 3))

g_pred <- ggplot(ge_test_obs, aes(x = goodman_filtered_expression_uniqueCount, y = ge_pred_obs)) +
  geom_point(colour = "#999999") +
  geom_smooth(method = lm, color = "black") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Observed Pleiotropy", 
       y = "Predicted Pleiotropy", 
       title = "GAP Expression Observed") +
  annotate(geom="text", label = test2, x = r2_text_x, y = 4500, size = annotate_size, parse = TRUE, hjust = 0) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  ylim(0,5500) +
  xlim(0,5500) +
  labs(tag = "D")

# goodman expression Permuted------------------------------------------------------
ge_pred_permuted <- predict(object = goodman_expr_perm, data = ge_test_permuted)
ge_test_permuted$ge_pred_permuted <- ge_pred_permuted$predictions

test2 <- paste("~R^2==~", 
               round(cor(ge_test_permuted$pleiotropy_gpe, ge_test_permuted$ge_pred_permuted)^2, digits = 3))

h_pred <- ggplot(ge_test_permuted, aes(x = pleiotropy_gpe, y = ge_pred_permuted)) +
  geom_point(colour = "#999999") +
  geom_smooth(method = lm, color = "black") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Permuted Pleiotropy", 
       y = "Predicted Pleiotropy", 
       title = "GAP Expression Permuted") +
  annotate(geom="text", label = test2, x = r2_text_x, y = 4500, size = annotate_size, parse = TRUE, hjust = 0) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  ylim(0,5500) +
  xlim(0,5500) +
  labs(tag = "D")

# test combined plot-------------------------------------------------------------
# ggpubr::ggarrange(a_pred,b_pred,c_pred,d_pred,
#                   e_pred,f_pred,g_pred,h_pred,
#                   nrow = 4, ncol = 2,
#                   common.legend = TRUE, legend = "bottom", align = "v")


# Export to file-----------------------------------------------------------------
# ggsave("~/git_projects/pleiotropy/images/model1_random_forest_prediction_accuracy_observed_vs_permuted.png",
#        plot = ggpubr::ggarrange(a_pred,b_pred,c_pred,d_pred, e_pred,f_pred,g_pred,h_pred,
#                                 nrow = 4, ncol = 2, 
#                                 common.legend = TRUE, legend = "bottom", align = "v"),
#        width = 11,
#        height = 10, 
#        units = "in",
#        dpi = "retina")

# just observed
ggsave("~/git_projects/pleiotropy/images/model1_random_forest_prediction_accuracy_observed.png",
       plot = ggpubr::ggarrange(a_pred,c_pred,e_pred,g_pred,
                                nrow = 2, ncol = 2, 
                                common.legend = TRUE, legend = "bottom", align = "v"),
       width = 11.5,
       height = 8, 
       units = "in",
       dpi = "retina")

# just permuted
ggsave("~/git_projects/pleiotropy/images/model1_random_forest_prediction_accuracy_permuted.png",
       plot = ggpubr::ggarrange(b_pred,d_pred,f_pred,h_pred,
                                nrow = 2, ncol = 2, 
                                common.legend = TRUE, legend = "bottom", align = "v"),
       width = 11.5,
       height = 8, 
       units = "in",
       dpi = "retina")

