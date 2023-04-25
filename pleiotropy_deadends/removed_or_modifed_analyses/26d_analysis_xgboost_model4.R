# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch and Emily Yi
# Contact... mbb262@cornell.edu
# Date...... 2022-03-23
# Updated... 2022-05-23
#
# Description: Model 4 (XGBooost without GO)
# Using pleiotropy scores from number of unique traits
# combined with many biological and noise responses in a XGBoost model
# 
# Running 8 models, 4: observed data, 4: permuted data
# Each of the 4 models comes from the 4 population-trait categories
#
# Model:
# pleiotropy score ~ ATACseq count + max RNA + max protein + avg. GERP
#                   + interval type + average LD + interval size 
#                   + input SNP count + adjusted expression pleiotropy
# ---------------------------------------------------------------

# Load in packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(xgboost)
library(testthat)

# set seed
set.seed(2022)

# read in data aggregatied by aggregate_interval.data.R
# data_dir <- "/Volumes/merrittData1/pleiotropy/interval_data"
data_dir <- "~/git_projects/pleiotropy/data/"
all_interval_data <- data.table::fread(paste0(data_dir, "/all_interval_data_unique_filter.csv"))

# Replace NAs with -1 for GERP columns
all_interval_data$np_noperm_mean_gerp[is.na(all_interval_data$np_noperm_mean_gerp)] <- -1
all_interval_data$gp_noperm_mean_gerp[is.na(all_interval_data$gp_noperm_mean_gerp)] <- -1
all_interval_data$gm_noperm_mean_gerp[is.na(all_interval_data$gm_noperm_mean_gerp)] <- -1
all_interval_data$ge_noperm_mean_gerp[is.na(all_interval_data$ge_noperm_mean_gerp)] <- -1


# -----------------------------------------------------------------------------
#           Format data
# -----------------------------------------------------------------------------

# Melt data for permuted columns ----------------------------------------------

all_interval_data_melt_permuted <- all_interval_data %>% 
  select("rr_id", 
         "range_type", 
         "seqid",
         "interval_size",
         "average_r2_ld_nam",
         "average_r2_ld_282",
         "snp_count_nam_all",
         "snp_count_goodman_all",
         "np_noperm_mean_gerp",
         "gp_noperm_mean_gerp",
         "gm_noperm_mean_gerp",
         "ge_noperm_mean_gerp",
         "max_protein_expression_23tissues",
         "max_rna_expression_23tissues",
         "atac_seq_peak_count", "residual_goodman_filtered_expression_uniqueCount",
         "npp_unique_perm_1","npp_unique_perm_2","npp_unique_perm_3","npp_unique_perm_4","npp_unique_perm_5",
         "npp_unique_perm_6","npp_unique_perm_7","npp_unique_perm_8","npp_unique_perm_9","npp_unique_perm_10",
         "gpp_unique_perm_1","gpp_unique_perm_2","gpp_unique_perm_3","gpp_unique_perm_4","gpp_unique_perm_5",
         "gpp_unique_perm_6","gpp_unique_perm_7","gpp_unique_perm_8","gpp_unique_perm_9","gpp_unique_perm_10",
         "gpm_perm_1","gpm_perm_2","gpm_perm_3","gpm_perm_4","gpm_perm_5",
         "gpm_perm_6","gpm_perm_7","gpm_perm_8","gpm_perm_9","gpm_perm_10",
         "gpe_perm_1","gpe_perm_2","gpe_perm_3","gpe_perm_4","gpe_perm_5")


# Melt individual datasets, then merge 
nam_field_perm <- tidyr::pivot_longer(all_interval_data_melt_permuted,
                                      cols = npp_unique_perm_1:npp_unique_perm_10,
                                      names_to = "pop_trait_perm_no_perm_level",
                                      values_to = "pleiotropy_npp")
goodman_field_perm <- tidyr::pivot_longer(all_interval_data_melt_permuted,
                                          cols = gpp_unique_perm_1:gpp_unique_perm_10,
                                          names_to = "pop_trait_perm_no_perm_level",
                                          values_to = "pleiotropy_gpp")
goodman_mass_perm <- tidyr::pivot_longer(all_interval_data_melt_permuted,
                                         cols = gpm_perm_1:gpm_perm_10,
                                         names_to = "pop_trait_perm_no_perm_level",
                                         values_to = "pleiotropy_gpm")

# Add a column with information about which permutation
nam_field_perm$perm_number <- gsub(".*_", "", nam_field_perm$pop_trait_perm_no_perm_level)
goodman_field_perm$perm_number <- gsub(".*_", "", goodman_field_perm$pop_trait_perm_no_perm_level)
goodman_mass_perm$perm_number <- gsub(".*_", "", goodman_mass_perm$pop_trait_perm_no_perm_level)

# Merge
permuted_merged <- Reduce(function(x, y) merge(x, y, by = c("rr_id", "perm_number", "range_type")),
                          list(nam_field_perm[,c(1:16,43:44)],
                               goodman_field_perm[,c(1:2,43:44)],
                               goodman_mass_perm[,c(1:2,43:44)]))
dim(permuted_merged) # should be: 754,900     20
colnames(permuted_merged)

# Keep this one separate because it only has 5 permutations
goodman_expression_perm <- tidyr::pivot_longer(all_interval_data_melt_permuted,
                                          cols = gpe_perm_1:gpe_perm_5,
                                          names_to = "pop_trait_perm_no_perm_level",
                                          values_to = "pleiotropy_gpe")

# Remove variables from space
rm(nam_field_perm)
rm(goodman_field_perm)
rm(goodman_mass_perm)


# ------------------------------------------------------------------------------
# Pipeline separated by population
# 1. Partition data into Observed Data and Permuted sets
#   and separate population-trait categories
# 2. Split training and testing data
#   FILTER BASED ON CHROMOSOME, 1:9 TRAIN, 10 TEST
# 3. Run xgboost random forest
# ------------------------------------------------------------------------------

# -----------------------------------------------------------
# NAM Observed Data
# -----------------------------------------------------------

# Partition data
np_obs <- all_interval_data %>% select("seqid",
                                        "average_r2_ld_nam",
                                        "np_noperm_mean_gerp",
                                        "max_protein_expression_23tissues",
                                        "max_rna_expression_23tissues",
                                        "atac_seq_peak_count",
                                        "snp_count_nam_all",
                                        "range_type",
                                        "interval_size",
                                        "nam_filtered_all_true_uniqueCount", 
                                        "residual_goodman_filtered_expression_uniqueCount")
dim(np_obs) # should be: [1] 75490    11
np_obs$range_type <- as.factor(np_obs$range_type)

# Split training and testing data
np_train_obs <- np_obs %>% filter(seqid != 10)
np_test_obs <- np_obs %>% filter(seqid == 10)

# prepare data for random forest - remove metadata, response
npf_train_data_mat <- np_train_obs[,-c("seqid","nam_filtered_all_true_uniqueCount")] %>% data.matrix()
npf_train_label_mat <- np_train_obs[,"nam_filtered_all_true_uniqueCount"] %>% data.matrix()
npf_test_data_mat <- np_test_obs[,-c("seqid","nam_filtered_all_true_uniqueCount")] %>% data.matrix()
npf_test_label_mat <- np_test_obs[,"nam_filtered_all_true_uniqueCount"] %>% data.matrix()
npf_dtrain <- xgboost::xgb.DMatrix(data = npf_train_data_mat, label=npf_train_label_mat)
npf_dtest <- xgboost::xgb.DMatrix(data = npf_test_data_mat, label=npf_test_label_mat)
npf_watchlist <- list(train=npf_dtrain, test=npf_dtest)

# run boosting
# compare hyperparameters
#============================
# sink("boosting_hyperparam.txt", append=FALSE, split=TRUE)
# for (et in 0.1*(1:10)) {
#   for (depth in 0:8) {
#       print(paste0("eta = ", et, " max.depth = ", depth))
#       xgb.train(data=dtrain, max.depth=depth, eta=et, nthread = 2, nrounds=10, watchlist=watchlist)
#   }
# }
# sink() # output back to terminal
#=============================
npf_boost <- xgb.train(data=npf_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=npf_watchlist)
npf_pred <- predict(npf_boost, npf_dtest)
npf_importance <- xgb.importance(model = npf_boost)

# clean up variables
rm(np_obs)
rm(npf_train_data_mat)
rm(npf_train_label_mat)
rm(npf_test_data_mat)
rm(np_train_obs)
rm(np_test_obs)
rm(npf_boost)
rm(npf_watchlist)
rm(npf_dtrain)
rm(npf_dtest)


# ------------------------------------------------------------
# NAM Permuted
# ------------------------------------------------------------

# prepare data for random forest - remove metadata, response, other populations
np_permuted <- permuted_merged %>% select("seqid",
                                          "range_type", 
                                          "interval_size",
                                          "average_r2_ld_nam",
                                          "snp_count_nam_all",
                                          "np_noperm_mean_gerp",
                                          "max_protein_expression_23tissues",
                                          "max_rna_expression_23tissues",
                                          "atac_seq_peak_count", 
                                          "residual_goodman_filtered_expression_uniqueCount",
                                          "pleiotropy_npp")
dim(np_permuted) # should be: [1] 754900     11
np_permuted$range_type <- as.factor(np_permuted$range_type)

np_train_permuted <- np_permuted %>% filter(seqid != 10)
np_test_permuted <- np_permuted %>% filter(seqid == 10)

# Remove seqid, response
npr_train_data_mat <- np_train_permuted %>% select(-c("seqid", "pleiotropy_npp")) %>% data.matrix()
npr_train_label_mat <- np_train_permuted[,"pleiotropy_npp"] %>% data.matrix()
npr_test_data_mat <- np_test_permuted %>% select(-c("seqid", "pleiotropy_npp")) %>% data.matrix()
npr_test_label_mat <- np_test_permuted[,"pleiotropy_npp"] %>% data.matrix()
npr_dtrain <- xgboost::xgb.DMatrix(data = npr_train_data_mat, label=npr_train_label_mat)
npr_dtest <- xgboost::xgb.DMatrix(data = npr_test_data_mat, label=npr_test_label_mat)
npr_watchlist <- list(train=npr_dtrain, test=npr_dtest)
npr_boost <- xgboost::xgb.train(data=npr_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=npr_watchlist)
npr_pred <- predict(npr_boost, npr_dtest)
npr_importance <- xgboost::xgb.importance(model = npr_boost)


# clean up variables
rm(np_permuted)
rm(npr_train_data_mat)
rm(npr_train_label_mat)
rm(npr_test_data_mat)
rm(np_train_permuted)
rm(np_test_permuted)
rm(npr_boost)
rm(npr_watchlist)
rm(npr_dtrain)
rm(npr_dtest)


# -----------------------------------------------------------
# goodman Field Observed Data
# -----------------------------------------------------------

gp_obs <- all_interval_data %>% select("seqid",
                                        "average_r2_ld_282",
                                        "gp_noperm_mean_gerp",
                                        "max_protein_expression_23tissues",
                                        "max_rna_expression_23tissues",
                                        "atac_seq_peak_count",
                                        "snp_count_goodman_all",
                                        "range_type",
                                        "interval_size",
                                        "goodman_filtered_physiological_true_uniqueCount", 
                                        "residual_goodman_filtered_expression_uniqueCount")
dim(gp_obs) # should be: [1] 75490    11
gp_obs$range_type <- as.factor(gp_obs$range_type)

gp_train_obs <- gp_obs %>% filter(seqid != 10)
gp_test_obs <- gp_obs %>% filter(seqid == 10)

# prepare data
gpf_train_data_mat <- gp_train_obs[,-c("seqid",
                                       "goodman_filtered_physiological_true_uniqueCount")] %>% data.matrix()
gpf_train_label_mat <- gp_train_obs[,"goodman_filtered_physiological_true_uniqueCount"] %>% data.matrix()
gpf_test_data_mat <- gp_test_obs[,-c("seqid",
                                     "goodman_filtered_physiological_true_uniqueCount")] %>% data.matrix()
gpf_test_label_mat <- gp_test_obs[,"goodman_filtered_physiological_true_uniqueCount"] %>% data.matrix()
gpf_dtrain <- xgboost::xgb.DMatrix(data = gpf_train_data_mat, label=gpf_train_label_mat)
gpf_dtest <- xgboost::xgb.DMatrix(data = gpf_test_data_mat, label=gpf_test_label_mat)
gpf_watchlist <- list(train=gpf_dtrain, test=gpf_dtest)

# run boosting
gpf_boost <- xgboost::xgb.train(data=gpf_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=gpf_watchlist)
gpf_pred <- predict(gpf_boost, gpf_dtest)
gpf_importance <- xgboost::xgb.importance(model = gpf_boost)

# clean up variables
rm(gp_obs)
rm(gpf_train_data_mat)
rm(gpf_train_label_mat)
rm(gpf_test_data_mat)
rm(gp_train_obs)
rm(gp_test_obs)
rm(gpf_boost)
rm(gpf_watchlist)
rm(gpf_dtrain)
rm(gpf_dtest)


# -----------------------------------------------------------
# Goodman Field Permuted
# -----------------------------------------------------------

gp_permuted <- permuted_merged %>% select("seqid",
                                          "range_type", 
                                          "interval_size",
                                          "average_r2_ld_282",
                                          "snp_count_goodman_all",
                                          "gp_noperm_mean_gerp",
                                          "max_protein_expression_23tissues",
                                          "max_rna_expression_23tissues",
                                          "atac_seq_peak_count",
                                          "pleiotropy_gpp", 
                                          "residual_goodman_filtered_expression_uniqueCount")
dim(gp_permuted) # should be: [1] 754900     11
gp_permuted$range_type <- as.factor(gp_permuted$range_type)

gp_train_permuted <- gp_permuted %>% filter(seqid != 10)
gp_test_permuted <- gp_permuted %>% filter(seqid == 10)

# prepare data
gpr_train_data_mat <- gp_train_permuted %>% select(-c("seqid", "pleiotropy_gpp")) %>% data.matrix()
gpr_train_label_mat <- gp_train_permuted[,"pleiotropy_gpp"] %>% data.matrix()
gpr_test_data_mat <- gp_test_permuted %>% select(-c("seqid", "pleiotropy_gpp")) %>% data.matrix()
gpr_test_label_mat <- gp_test_permuted[,"pleiotropy_gpp"] %>% data.matrix()
gpr_dtrain <- xgboost::xgb.DMatrix(data = gpr_train_data_mat, label=gpr_train_label_mat)
gpr_dtest <- xgboost::xgb.DMatrix(data = gpr_test_data_mat, label=gpr_test_label_mat)
gpr_watchlist <- list(train=gpr_dtrain, test=gpr_dtest)

# run boosting
gpr_boost <- xgboost::xgb.train(data=gpr_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=gpr_watchlist)
gpr_pred <- predict(gpr_boost, gpr_dtest)
gpr_importance <- xgboost::xgb.importance(model = gpr_boost)

# clean up variables
rm(gp_permuted)
rm(gpr_train_data_mat)
rm(gpr_train_label_mat)
rm(gpr_test_data_mat)
rm(gp_train_permuted)
rm(gp_test_permuted)
rm(gpr_boost)
rm(gpr_watchlist)
rm(gpr_dtrain)
rm(gpr_dtest)


# -----------------------------------------------------------
# Goodman Mass Features Observed Data
# -----------------------------------------------------------

gm_obs <- all_interval_data %>% select("seqid",
                                        "average_r2_ld_282",
                                        "gm_noperm_mean_gerp",
                                        "max_protein_expression_23tissues",
                                        "max_rna_expression_23tissues",
                                        "atac_seq_peak_count",
                                        "snp_count_goodman_all",
                                        "range_type",
                                        "interval_size",
                                        "goodman_filtered_metabolite_true_uniqueCount", 
                                        "residual_goodman_filtered_expression_uniqueCount")
dim(gm_obs) # should be: [1] 75490    11
gm_obs$range_type <- as.factor(gm_obs$range_type)

gm_train_obs <- gm_obs %>% filter(seqid != 10)
gm_test_obs <- gm_obs %>% filter(seqid == 10)

# prepare data
gmf_train_data_mat <- gm_train_obs[,-c("seqid", "goodman_filtered_metabolite_true_uniqueCount")] %>% data.matrix()
gmf_train_label_mat <- gm_train_obs[,"goodman_filtered_metabolite_true_uniqueCount"] %>% data.matrix()
gmf_test_data_mat <- gm_test_obs[,-c("seqid", "goodman_filtered_metabolite_true_uniqueCount")] %>% data.matrix()
gmf_test_label_mat <- gm_test_obs[,"goodman_filtered_metabolite_true_uniqueCount"] %>% data.matrix()
gmf_dtrain <- xgboost::xgb.DMatrix(data = gmf_train_data_mat, label=gmf_train_label_mat)
gmf_dtest <- xgboost::xgb.DMatrix(data = gmf_test_data_mat, label=gmf_test_label_mat)
gmf_watchlist <- list(train=gmf_dtrain, test=gmf_dtest)

# run boosting
# best at max.depth=7
gmf_boost <- xgboost::xgb.train(data=gmf_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=gmf_watchlist)
gmf_pred <- predict(gmf_boost, gmf_dtest)
gmf_importance <- xgboost::xgb.importance(model = gmf_boost)

# clean up variables
rm(gm_obs)
rm(gmf_train_data_mat)
rm(gmf_train_label_mat)
rm(gmf_test_data_mat)
rm(gm_train_obs)
rm(gm_test_obs)
rm(gmf_boost)
rm(gmf_watchlist)
rm(gmf_dtrain)
rm(gmf_dtest)


# -----------------------------------------------------------
# Goodman Mass Features Permuted
# -----------------------------------------------------------

gm_permuted <- permuted_merged %>% select("seqid",
                                          "range_type", 
                                          "interval_size",
                                          "average_r2_ld_282",
                                          "snp_count_goodman_all",
                                          "gm_noperm_mean_gerp",
                                          "max_protein_expression_23tissues",
                                          "max_rna_expression_23tissues",
                                          "atac_seq_peak_count",
                                          "pleiotropy_gpm", 
                                          "residual_goodman_filtered_expression_uniqueCount")
dim(gm_permuted) # should be: [1] 754900     11
gm_permuted$range_type <- as.factor(gm_permuted$range_type)

gm_train_permuted <- gm_permuted %>% filter(seqid != 10)
gm_test_permuted <- gm_permuted %>% filter(seqid == 10)

# prepare data
gmr_train_data_mat <- gm_train_permuted %>% select(-c("seqid", "pleiotropy_gpm")) %>% data.matrix()
gmr_train_label_mat <- gm_train_permuted[,"pleiotropy_gpm"] %>% data.matrix()
gmr_test_data_mat <- gm_test_permuted %>% select(-c("seqid", "pleiotropy_gpm")) %>% data.matrix()
gmr_test_label_mat <- gm_test_permuted[,"pleiotropy_gpm"] %>% data.matrix()
gmr_dtrain <- xgboost::xgb.DMatrix(data = gmr_train_data_mat, label=gmr_train_label_mat)
gmr_dtest <- xgboost::xgb.DMatrix(data = gmr_test_data_mat, label=gmr_test_label_mat)
gmr_watchlist <- list(train=gmr_dtrain, test=gmr_dtest)

# run boosting
gmr_boost <- xgboost::xgb.train(data=gmr_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=gmr_watchlist)
gmr_pred <- predict(gmr_boost, gmr_dtest)
gmr_importance <- xgboost::xgb.importance(model = gmr_boost)


# clean up variables
rm(gm_permuted)
rm(gmr_train_data_mat)
rm(gmr_train_label_mat)
rm(gmr_test_data_mat)
rm(gm_train_permuted)
rm(gm_test_permuted)
rm(gmr_boost)
rm(gmr_watchlist)
rm(gmr_dtrain)
rm(gmr_dtest)


# ------------------------------------------------------------
# goodman expression Observed Data
# ------------------------------------------------------------

ge_obs <- all_interval_data %>% select("seqid",
                                        "average_r2_ld_282",
                                        "ge_noperm_mean_gerp",
                                        "max_protein_expression_23tissues",
                                        "max_rna_expression_23tissues",
                                        "atac_seq_peak_count",
                                        "snp_count_goodman_all",
                                        "range_type",
                                        "interval_size",
                                        "goodman_filtered_expression_uniqueCount")
dim(ge_obs) # should be: [1] 75490    11
ge_obs$range_type <- as.factor(ge_obs$range_type)

ge_train_obs <- ge_obs %>% filter(seqid != 10)
ge_test_obs <- ge_obs %>% filter(seqid == 10)

# prepare data
gef_train_data_mat <- ge_train_obs[,-c("seqid","goodman_filtered_expression_uniqueCount")] %>% data.matrix()
gef_train_label_mat <- ge_train_obs[,"goodman_filtered_expression_uniqueCount"] %>% data.matrix()
gef_test_data_mat <- ge_test_obs[,-c("seqid","goodman_filtered_expression_uniqueCount")] %>% data.matrix()
gef_test_label_mat <- ge_test_obs[,"goodman_filtered_expression_uniqueCount"] %>% data.matrix()
gef_dtrain <- xgboost::xgb.DMatrix(data = gef_train_data_mat, label=gef_train_label_mat)
gef_dtest <- xgboost::xgb.DMatrix(data = gef_test_data_mat, label=gef_test_label_mat)
gef_watchlist <- list(train=gef_dtrain, test=gef_dtest)

# run boosting
gef_boost <- xgboost::xgb.train(data=gef_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=gef_watchlist)
gef_pred <- predict(gef_boost, gef_dtest)
gef_importance <- xgboost::xgb.importance(model = gef_boost)

# clean up variables
rm(ge_obs)
rm(gef_train_data_mat)
rm(gef_train_label_mat)
rm(gef_test_data_mat)
rm(ge_train_obs)
rm(ge_test_obs)
rm(gef_boost)
rm(gef_watchlist)
rm(gef_dtrain)
rm(gef_dtest)


# -----------------------------------------------------------
# Goodman Expression Permuted
# -----------------------------------------------------------

ge_permuted <- goodman_expression_perm %>% select("seqid",
                                                  "range_type", 
                                                  "interval_size",
                                                  "average_r2_ld_282",
                                                  "snp_count_goodman_all",
                                                  "ge_noperm_mean_gerp",
                                                  "max_protein_expression_23tissues",
                                                  "max_rna_expression_23tissues",
                                                  "atac_seq_peak_count",
                                                  "pleiotropy_gpe")
dim(ge_permuted) # should be: [1] 377450     11
ge_permuted$range_type <- as.factor(ge_permuted$range_type)

ge_train_permuted <- ge_permuted %>% filter(seqid != 10)
ge_test_permuted <- ge_permuted %>% filter(seqid == 10)

# prepare data
ger_train_data_mat <- ge_train_permuted %>% select(-c("seqid", "pleiotropy_gpe")) %>% data.matrix()
ger_train_label_mat <- ge_train_permuted[,"pleiotropy_gpe"] %>% data.matrix()
ger_test_data_mat <- ge_test_permuted %>% select(-c("seqid", "pleiotropy_gpe")) %>% data.matrix()
ger_test_label_mat <- ge_test_permuted[,"pleiotropy_gpe"] %>% data.matrix()
ger_dtrain <- xgb.DMatrix(data = ger_train_data_mat, label=ger_train_label_mat)
ger_dtest <- xgb.DMatrix(data = ger_test_data_mat, label=ger_test_label_mat)
ger_watchlist <- list(train=ger_dtrain, test=ger_dtest)

# run boosting
ger_boost <- xgb.train(data=ger_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=ger_watchlist)
ger_pred <- predict(ger_boost, ger_dtest)
ger_importance <- xgb.importance(model = ger_boost)

# clean up variables
rm(ge_permuted)
rm(ger_train_data_mat)
rm(ger_train_label_mat)
rm(ger_test_data_mat)
rm(ge_train_permuted)
rm(ge_test_permuted)
rm(ger_boost)
rm(ger_watchlist)
rm(ger_dtrain)
rm(ger_dtest)


# --------------------------------------------------------------------------------
# Make a table summarizing results
# --------------------------------------------------------------------------------

# Returns MSE and R^2
fit.stats <- function(pred, labels) {
  c(mean((pred - labels)^2), cor(pred, labels)^2) %>% return
}

# non-adjusted data, observed and permuted
models <- c("NAM Field Observed Data", 
            "NAM Field Permuted Data",
            "GAP Field Observed Data",
            "GAP Field Permuted Data",
            "GAP Mass Features Observed Data",
            "GAP Mass Features Permuted Data",
            "GAP Expression Observed Data",
            "GAP Expression Permuted Data")
temp <- data.frame(matrix(ncol = 3, nrow = 8))
colnames(temp) <- c("model", "MSE", "r_squared")
temp$model <- models
temp[1,2:3] <- fit.stats(npf_pred, npf_test_label_mat)
temp[2,2:3] <- fit.stats(npr_pred, npr_test_label_mat)
temp[3,2:3] <- fit.stats(gpf_pred, gpf_test_label_mat)
temp[4,2:3] <- fit.stats(gpr_pred, gpr_test_label_mat)
temp[5,2:3] <- fit.stats(gmf_pred, gmf_test_label_mat)
temp[6,2:3] <- fit.stats(gmr_pred, gmr_test_label_mat)
temp[7,2:3] <- fit.stats(gef_pred, gef_test_label_mat)
temp[8,2:3] <- fit.stats(ger_pred, ger_test_label_mat)

# Export
data.table::fwrite(temp, "~/git_projects/pleiotropy/data/model4_xgboost_prediction_accuracy.csv")


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

# Gather nice names for plots
pretty_names <- data.table::fread("~/git_projects/pleiotropy/data/pretty_names_for_rf.csv")

prepare_importance <- function(old) {
  new <- as.data.frame(old)
  names(new)[1:2] <- c("Variable", "Importance")
  new$Importance <- new$Importance/max(new$Importance)
  new <- merge(x = new, y = pretty_names, 
               by.x = "Variable", by.y = "old", all.x = TRUE)
  new <- new %>%
    select(-c("Variable", "Cover", "Frequency"))
  names(new)[2] <- "Variable"
  return(new)
}


# NAM Phys Observed Data -----------------------------------------------
npf_importance <- prepare_importance(npf_importance)
a <- npf_importance %>% 
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
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9","#E69F00")) +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
  labs(tag = "A")

# nam phys Permuted -------------------------------------------------
npr_importance <- prepare_importance(npr_importance)
b <- npr_importance %>% 
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
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9","#E69F00")) +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
  labs(tag = "B")

# goodman Field Observed Data ----------------------------------
gpf_importance <- prepare_importance(gpf_importance)
c <- gpf_importance %>% 
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
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9","#E69F00")) +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
  labs(tag = "C")

# goodman Field Permuted ------------------------------------
gpr_importance <- prepare_importance(gpr_importance)
d <- gpr_importance %>% 
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
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9","#E69F00")) +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
  labs(tag = "D")

# goodman Mass Features Observed Data -------------------------------------
gmf_importance<- prepare_importance(gmf_importance)
e <- gmf_importance %>% 
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
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9","#E69F00")) +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
  labs(tag = "E")

# goodman Mass Features Permuted --------------------------------------
gmr_importance <- prepare_importance(gmr_importance)
f <- gmr_importance %>% 
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
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9","#E69F00")) +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
  labs(tag = "F")

# goodman expression Observed Data ------------------------------------
gef_importance <- prepare_importance(gef_importance)
g <- gef_importance %>% 
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
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=c("#999999", "#56B4E9","#E69F00")) +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
  labs(tag = "G")

# goodman expression Permuted --------------------------------------
ger_importance <- prepare_importance(ger_importance)
h <- ger_importance %>% 
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
        plot.tag=element_text(size=tag_size),
        legend.position="none") +
  scale_colour_manual(values=c("#999999", "#56B4E9","#E69F00")) +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00")) +
  labs(tag = "H")


# test combined plot ----------------------------------------------
# ggpubr::ggarrange(a,b,c,d,
#                   e,f,g,h,
#                   nrow = 4, ncol = 2,
#                   common.legend = TRUE, legend = "bottom", align = "v")

# Export to file---------------------------------------------------
# Observed data
ggsave("~/git_projects/pleiotropy/images/model4_random_forest_importance_observed.png",
       plot = ggpubr::ggarrange(a,c,e,g, nrow = 2, ncol = 2, 
                                common.legend = TRUE, legend = "bottom", align = "v"),
       width = 17,
       height = 13, 
       units = "in",
       dpi = "retina")

# Permuted data
ggsave("~/git_projects/pleiotropy/images/model4_random_forest_importance_permuted.png",
       plot = ggpubr::ggarrange(b,d,f,h,nrow = 2, ncol = 2, 
                                common.legend = TRUE, legend = "bottom", align = "v"),
       width = 17,
       height = 13, 
       units = "in",
       dpi = "retina")

# Combine observed and permuted in one plot
ggsave("~/git_projects/pleiotropy/images/model4_random_forest_importance_observed_vs_permuted.png",
       plot = ggpubr::ggarrange(a,b,c,d, e,f,g,h,
                                nrow = 4, ncol = 2, 
                                common.legend = TRUE, legend = "bottom", align = "v"),
       width = 16,
       height = 15, 
       units = "in",
       dpi = "retina")


# --------------------------------------------------------------------------------
# Plot & combine all model performance plots 
# --------------------------------------------------------------------------------

# nam phys Observed Data --------------------------------------------------------------
np_test_obs = data.frame(nam_filtered_all_true_uniqueCount = npf_test_label_mat, np_pred_obs = npf_pred) 
test2 <- paste("~R^2==~", round(temp$r_squared[1], digits = 3))

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
        plot.tag=element_text(size=tag_size)) +
  ylim(-1,65) +
  xlim(0,65) +
  labs(tag = "A")


# nam phys Permuted-----------------------------------------------------------------
np_test_permuted = data.frame(pleiotropy_npp = npr_test_label_mat, np_pred_permuted = npr_pred)
test2 <- paste("~R^2==~", round(temp$r_squared[2], digits = 3))

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
        plot.tag=element_text(size=tag_size)) +
  ylim(-1,65) +
  xlim(0,65) +
  labs(tag = "B")

# goodman Field Observed Data--------------------------------------------------
gp_test_obs = data.frame(goodman_filtered_physiological_true_uniqueCount = gpf_test_label_mat, gp_pred_obs = gpf_pred)
test2 <- paste("~R^2==~", round(temp$r_squared[3], digits = 3))

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
        plot.tag=element_text(size=tag_size)) +
  ylim(-1,65) +
  xlim(0,65) +
  labs(tag = "C")

# goodman Field Permuted---------------------------------------------------
gp_test_permuted = data.frame(pleiotropy_gpp = gpr_test_label_mat, gp_pred_permuted = gpr_pred)
test2 <- paste("~R^2==~", round(temp$r_squared[4], digits = 3))

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
        plot.tag=element_text(size=tag_size)) +
  ylim(-1,65) +
  xlim(0,65) +
  labs(tag = "D")

# goodman Mass Features Observed Data----------------------------------------------------
gm_test_obs = data.frame(goodman_filtered_metabolite_true_uniqueCount = gmf_test_label_mat, gm_pred_obs = gmf_pred)
test2 <- paste("~R^2==~", round(temp$r_squared[5], digits = 3))

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
        plot.tag=element_text(size=tag_size)) +
  ylim(-32,1850) +
  xlim(0,1850) +
  labs(tag = "E")

# goodman Mass Features Permuted-----------------------------------------------------
gm_test_permuted = data.frame(pleiotropy_gpm = gmr_test_label_mat, gm_pred_permuted = gmr_pred)
test2 <- paste("~R^2==~", round(temp$r_squared[6], digits = 3))

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
        plot.tag=element_text(size=tag_size)) +
  ylim(-32,1850) +
  xlim(0,1850) +
  labs(tag = "F")

# goodman expression Observed Data---------------------------------------------------
ge_test_obs = data.frame(goodman_filtered_expression_uniqueCount = gef_test_label_mat, ge_pred_obs = gef_pred)
test2 <- paste("~R^2==~", round(temp$r_squared[7], digits = 3))

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
        plot.tag=element_text(size=tag_size)) +
  ylim(-18,5500) +
  xlim(0,5500) +
  labs(tag = "G")

# goodman expression Permuted------------------------------------------------------
ge_test_permuted = data.frame(pleiotropy_gpe = ger_test_label_mat, ge_pred_permuted = ger_pred)
test2 <- paste("~R^2==~", round(temp$r_squared[8], digits = 3))

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
        plot.tag=element_text(size=tag_size)) +
  ylim(-18,5500) +
  xlim(0,5500) +
  labs(tag = "H")

# test combined plot-------------------------------------------------------------
ggpubr::ggarrange(a_pred,b_pred,c_pred,d_pred,
                  e_pred,f_pred,g_pred,h_pred,
                  nrow = 4, ncol = 2,
                  common.legend = TRUE, legend = "bottom", align = "v")


# Export to file-----------------------------------------------------------------
ggsave("~/git_projects/pleiotropy/images/model4_random_forest_prediction_accuracy_observed_vs_permuted.png",
       plot = ggpubr::ggarrange(a_pred,b_pred,c_pred,d_pred, e_pred,f_pred,g_pred,h_pred,
                                nrow = 4, ncol = 2, 
                                common.legend = TRUE, legend = "bottom", align = "v"),
       width = 11,
       height = 10, 
       units = "in",
       dpi = "retina")

