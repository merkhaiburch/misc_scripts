#XGBoost but for adjusted counts

# ---------------------------------------------------------------------------
# Repeat formatting for permuted adjusted (residual) counts
# ---------------------------------------------------------------------------

res_cols = (all_interval_data_melt_permuted %>% names)[1:15] %>%
  c((all_interval_data_melt_permuted %>% names)[16:50] %>% paste0("residual_", .))
all_interval_data_melt_res <- all_interval_data %>% select(all_of(res_cols))
rm(res_cols)

# Melt individual datasets, then merge 
nam_field_res <- tidyr::pivot_longer(all_interval_data_melt_res,
                                     cols = residual_npp_perm_1:residual_npp_perm_10,
                                     names_to = "pop_trait_perm_no_perm_level",
                                     values_to = "pleiotropy_npp")
goodman_field_res <- tidyr::pivot_longer(all_interval_data_melt_res,
                                         cols = residual_gpp_perm_1:residual_gpp_perm_10,
                                         names_to = "pop_trait_perm_no_perm_level",
                                         values_to = "pleiotropy_gpp")
goodman_mass_res <- tidyr::pivot_longer(all_interval_data_melt_res,
                                        cols = residual_gpm_perm_1:residual_gpm_perm_10,
                                        names_to = "pop_trait_perm_no_perm_level",
                                        values_to = "pleiotropy_gpm")

# Add a column with information about which permutation
nam_field_res$perm_number <- gsub(".*_", "", nam_field_res$pop_trait_perm_no_perm_level)
goodman_field_res$perm_number <- gsub(".*_", "", goodman_field_res$pop_trait_perm_no_perm_level)
goodman_mass_res$perm_number <- gsub(".*_", "", goodman_mass_res$pop_trait_perm_no_perm_level)

# Merge
res_merged <- Reduce(function(x, y) merge(x, y, by = c("rr_id", "perm_number")), 
                     list(nam_field_res[,c(1:15, 40:43)], 
                          goodman_field_res[,c(1, 40:43)], 
                          goodman_mass_res[,c(1, 40:43)]))
dim(res_merged) # should be: 754900     25

# Keep this one separate because it only has 5 permutations
goodman_expression_res <- tidyr::pivot_longer(all_interval_data_melt_res,
                                              cols = residual_gpe_perm_1:residual_gpe_perm_5,
                                              names_to = "pop_trait_perm_no_perm_level",
                                              values_to = "pleiotropy_gpe")
rm(nam_field_res)
rm(goodman_field_res)
rm(goodman_mass_res)


# --------------------------
# NAM adjusted counts
# --------------------------
# Partition data, the same as unadjusted except response is residual_...
np_res_filt <- all_interval_data %>% select("seqid",
                                            "average_r2_ld_nam",
                                            "np_noperm_mean_gerp",
                                            "max_protein_expression_23tissues",
                                            "max_rna_expression_23tissues",
                                            "atac_seq_peak_count",
                                            "snp_count_nam_all",
                                            "range_type",
                                            "interval_size",
                                            "residual_nam_filtered_all_uniqueCount")

np_res_filt$range_type <- as.factor(np_res_filt$range_type)

np_res_train_filt <- np_res_filt %>% filter(seqid != 10)
np_res_test_filt <- np_res_filt %>% filter(seqid == 10)

npf_res_train_data_mat <- np_res_train_filt[,-c("seqid",
                                                "residual_nam_filtered_all_uniqueCount"
)] %>% data.matrix()()
npf_res_train_label_mat <- np_res_train_filt[,"residual_nam_filtered_all_uniqueCount"] %>% data.matrix()
npf_res_test_data_mat <- np_res_test_filt[,-c("seqid",
                                              "residual_nam_filtered_all_uniqueCount"
)] %>% data.matrix()
npf_res_test_label_mat <- np_res_test_filt[,"residual_nam_filtered_all_uniqueCount"] %>% data.matrix()
npf_res_dtrain <- xgb.DMatrix(data = npf_res_train_data_mat, label=npf_res_train_label_mat)
npf_res_dtest <- xgb.DMatrix(data = npf_res_test_data_mat, label=npf_res_test_label_mat)
npf_res_watchlist <- list(train=npf_res_dtrain, test=npf_res_dtest)

npf_res_boost <- xgb.train(data=npf_res_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=npf_res_watchlist)
npf_res_pred <- predict(npf_res_boost, npf_res_dtest)
npf_res_importance <- xgb.importance(model = npf_res_boost)

# clean up variables
rm(np_res_filt)
rm(npf_res_train_data_mat)
rm(npf_res_train_label_mat)
rm(npf_res_test_data_mat)
rm(np_res_train_filt)
rm(np_res_test_filt)
rm(npf_res_boost)
rm(npf_res_watchlist)
rm(npf_res_dtrain)
rm(npf_res_dtest)

# ------------------------------
# NAM adjusted (residual) counts, permuted
# ------------------------------
# prepare data for random forest - remove metadata, response, other populations
np_res_permuted <- res_merged %>% select("seqid",
                                         "range_type", 
                                         "interval_size",
                                         "average_r2_ld_nam",
                                         "snp_count_nam_all",
                                         "np_noperm_mean_gerp",
                                         "max_protein_expression_23tissues",
                                         "max_rna_expression_23tissues",
                                         "atac_seq_peak_count",
                                         "pleiotropy_npp")
np_res_permuted$range_type <- as.factor(np_res_permuted$range_type)

np_res_train_permuted <- np_res_permuted %>% filter(seqid != 10)
np_res_test_permuted <- np_res_permuted %>% filter(seqid == 10)

# Remove seqid, response
npr_res_train_data_mat <- np_res_train_permuted %>%
  select(-c("seqid", "pleiotropy_npp")) %>% data.matrix()
npr_res_train_label_mat <- np_res_train_permuted[,"pleiotropy_npp"] %>% data.matrix()
npr_res_test_data_mat <- np_res_test_permuted %>%
  select(-c("seqid", "pleiotropy_npp")) %>% data.matrix()
npr_res_test_label_mat <- np_res_test_permuted[,"pleiotropy_npp"] %>% data.matrix()
npr_res_dtrain <- xgb.DMatrix(data = npr_res_train_data_mat, label=npr_res_train_label_mat)
npr_res_dtest <- xgb.DMatrix(data = npr_res_test_data_mat, label=npr_res_test_label_mat)
npr_res_watchlist <- list(train=npr_res_dtrain, test=npr_res_dtest)
npr_res_boost <- xgb.train(data=npr_res_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=npr_res_watchlist)
npr_res_pred <- predict(npr_res_boost, npr_res_dtest)
npr_res_importance <- xgb.importance(model = npr_res_boost)


# clean up variables
rm(np_res_permuted)
rm(npr_res_train_data_mat)
rm(npr_res_train_label_mat)
rm(npr_res_test_data_mat)
rm(np_res_train_permuted)
rm(np_res_test_permuted)
rm(npr_res_boost)
rm(npr_res_watchlist)
rm(npr_res_dtrain)
rm(npr_res_dtest)


# --------------------------
# Goodman Field (physiological) adjusted counts
# --------------------------
gp_res_filt <- all_interval_data %>% select("seqid",
                                            "average_r2_ld_282",
                                            "gp_noperm_mean_gerp",
                                            "max_protein_expression_23tissues",
                                            "max_rna_expression_23tissues",
                                            "atac_seq_peak_count",
                                            "snp_count_goodman_all",
                                            "range_type",
                                            "interval_size",
                                            "residual_goodman_filtered_physiological_uniqueCount")

gp_res_filt$range_type <- as.factor(gp_res_filt$range_type)

gp_res_train_filt <- gp_res_filt %>% filter(seqid != 10)
gp_res_test_filt <- gp_res_filt %>% filter(seqid == 10)

# prepare data
gpf_res_train_data_mat <- gp_res_train_filt[,-c("seqid",
                                                "residual_goodman_filtered_physiological_uniqueCount")] %>% data.matrix()
gpf_res_train_label_mat <- gp_res_train_filt[,"residual_goodman_filtered_physiological_uniqueCount"] %>% data.matrix()
gpf_res_test_data_mat <- gp_res_test_filt[,-c("seqid",
                                              "residual_goodman_filtered_physiological_uniqueCount")] %>% data.matrix()
gpf_res_test_label_mat <- gp_res_test_filt[,"residual_goodman_filtered_physiological_uniqueCount"] %>% data.matrix()
gpf_res_dtrain <- xgb.DMatrix(data = gpf_res_train_data_mat, label=gpf_res_train_label_mat)
gpf_res_dtest <- xgb.DMatrix(data = gpf_res_test_data_mat, label=gpf_res_test_label_mat)
gpf_res_watchlist <- list(train=gpf_res_dtrain, test=gpf_res_dtest)

# run boosting
gpf_res_boost <- xgb.train(data=gpf_res_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=gpf_res_watchlist)
gpf_res_pred <- predict(gpf_res_boost, gpf_res_dtest)
gpf_res_importance <- xgb.importance(model = gpf_res_boost)


# clean up variables
rm(gp_res_filt)
rm(gpf_res_train_data_mat)
rm(gpf_res_train_label_mat)
rm(gpf_res_test_data_mat)
rm(gp_res_train_filt)
rm(gp_res_test_filt)
rm(gpf_res_boost)
rm(gpf_res_watchlist)
rm(gpf_res_dtrain)
rm(gpf_res_dtest)

# ------------------------------
# Goodman Field adjusted (residual) counts, permuted
# ------------------------------
gp_res_permuted <- res_merged %>% select("seqid",
                                         "range_type", 
                                         "interval_size",
                                         "average_r2_ld_282",
                                         "snp_count_goodman_all",
                                         "gp_noperm_mean_gerp",
                                         "max_protein_expression_23tissues",
                                         "max_rna_expression_23tissues",
                                         "atac_seq_peak_count",
                                         "pleiotropy_gpp")


gp_res_permuted$range_type <- as.factor(gp_res_permuted$range_type)

gp_res_train_permuted <- gp_res_permuted %>% filter(seqid != 10)
gp_res_test_permuted <- gp_res_permuted %>% filter(seqid == 10)

# prepare data
gpr_res_train_data_mat <- gp_res_train_permuted %>%
  select(-c("seqid", "pleiotropy_gpp")) %>% data.matrix()
gpr_res_train_label_mat <- gp_res_train_permuted[,"pleiotropy_gpp"] %>% data.matrix()
gpr_res_test_data_mat <- gp_res_test_permuted %>%
  select(-c("seqid", "pleiotropy_gpp")) %>% data.matrix()
gpr_res_test_label_mat <- gp_res_test_permuted[,"pleiotropy_gpp"] %>% data.matrix()
gpr_res_dtrain <- xgb.DMatrix(data = gpr_res_train_data_mat, label=gpr_res_train_label_mat)
gpr_res_dtest <- xgb.DMatrix(data = gpr_res_test_data_mat, label=gpr_res_test_label_mat)
gpr_res_watchlist <- list(train=gpr_res_dtrain, test=gpr_res_dtest)

# run boosting
gpr_res_boost <- xgb.train(data=gpr_res_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=gpr_res_watchlist)
gpr_res_pred <- predict(gpr_res_boost, gpr_res_dtest)
gpr_res_importance <- xgb.importance(model = gpr_res_boost)


# clean up variables
rm(gp_res_permuted)
rm(gpr_res_train_data_mat)
rm(gpr_res_train_label_mat)
rm(gpr_res_test_data_mat)
rm(gp_res_train_permuted)
rm(gp_res_test_permuted)
rm(gpr_res_boost)
rm(gpr_res_watchlist)
rm(gpr_res_dtrain)
rm(gpr_res_dtest)


# --------------------------
# Goodman Mass Features adjusted counts
# --------------------------
gm_res_filt <- all_interval_data %>% select("seqid",
                                            "average_r2_ld_282",
                                            "gm_noperm_mean_gerp",
                                            "max_protein_expression_23tissues",
                                            "max_rna_expression_23tissues",
                                            "atac_seq_peak_count",
                                            "snp_count_goodman_all",
                                            "range_type",
                                            "interval_size",
                                            "residual_goodman_filtered_metabolite_uniqueCount")

gm_res_filt$range_type <- as.factor(gm_res_filt$range_type)

gm_res_train_filt <- gm_res_filt %>% filter(seqid != 10)
gm_res_test_filt <- gm_res_filt %>% filter(seqid == 10)

# prepare data
gmf_res_train_data_mat <- gm_res_train_filt[,-c("seqid",
                                                "residual_goodman_filtered_metabolite_uniqueCount")] %>% data.matrix()
gmf_res_train_label_mat <- gm_res_train_filt[,"residual_goodman_filtered_metabolite_uniqueCount"] %>% data.matrix()
gmf_res_test_data_mat <- gm_res_test_filt[,-c("seqid",
                                              "residual_goodman_filtered_metabolite_uniqueCount")] %>% data.matrix()
gmf_res_test_label_mat <- gm_res_test_filt[,"residual_goodman_filtered_metabolite_uniqueCount"] %>% data.matrix()
gmf_res_dtrain <- xgb.DMatrix(data = gmf_res_train_data_mat, label=gmf_res_train_label_mat)
gmf_res_dtest <- xgb.DMatrix(data = gmf_res_test_data_mat, label=gmf_res_test_label_mat)
gmf_res_watchlist <- list(train=gmf_res_dtrain, test=gmf_res_dtest)

# run boosting
gmf_res_boost <- xgb.train(data=gmf_res_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=gmf_res_watchlist)
gmf_res_pred <- predict(gmf_res_boost, gmf_res_dtest)
gmf_res_importance <- xgb.importance(model = gmf_res_boost)


# clean up variables
rm(gm_res_filt)
rm(gmf_res_train_data_mat)
rm(gmf_res_train_label_mat)
rm(gmf_res_test_data_mat)
rm(gm_res_train_filt)
rm(gm_res_test_filt)
rm(gmf_res_boost)
rm(gmf_res_watchlist)
rm(gmf_res_dtrain)
rm(gmf_res_dtest)

# ------------------------------
# Goodman Mass adjusted (residual) counts, permuted
# ------------------------------
gm_res_permuted <- res_merged %>% select("seqid",
                                         "range_type", 
                                         "interval_size",
                                         "average_r2_ld_282",
                                         "snp_count_goodman_all",
                                         "gm_noperm_mean_gerp",
                                         "max_protein_expression_23tissues",
                                         "max_rna_expression_23tissues",
                                         "atac_seq_peak_count",
                                         "pleiotropy_gpm")


gm_res_permuted$range_type <- as.factor(gm_res_permuted$range_type)

gm_res_train_permuted <- gm_res_permuted %>% filter(seqid != 10)
gm_res_test_permuted <- gm_res_permuted %>% filter(seqid == 10)

# prepare data
gmr_res_train_data_mat <- gm_res_train_permuted %>%
  select(-c("seqid", "pleiotropy_gpm")) %>% data.matrix()
gmr_res_train_label_mat <- gm_res_train_permuted[,"pleiotropy_gpm"] %>% data.matrix()
gmr_res_test_data_mat <- gm_res_test_permuted %>%
  select(-c("seqid", "pleiotropy_gpm")) %>% data.matrix()
gmr_res_test_label_mat <- gm_res_test_permuted[,"pleiotropy_gpm"] %>% data.matrix()
gmr_res_dtrain <- xgb.DMatrix(data = gmr_res_train_data_mat, label=gmr_res_train_label_mat)
gmr_res_dtest <- xgb.DMatrix(data = gmr_res_test_data_mat, label=gmr_res_test_label_mat)
gmr_res_watchlist <- list(train=gmr_res_dtrain, test=gmr_res_dtest)

# run boosting
gmr_res_boost <- xgb.train(data=gmr_res_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=gmr_res_watchlist)
gmr_res_pred <- predict(gmr_res_boost, gmr_res_dtest)
gmr_res_importance <- xgb.importance(model = gmr_res_boost)


# clean up variables
rm(gm_res_permuted)
rm(gmr_res_train_data_mat)
rm(gmr_res_train_label_mat)
rm(gmr_res_test_data_mat)
rm(gm_res_train_permuted)
rm(gm_res_test_permuted)
rm(gmr_res_boost)
rm(gmr_res_watchlist)
rm(gmr_res_dtrain)
rm(gmr_res_dtest)


# --------------------------
# Goodman Expression adjusted counts
# --------------------------
ge_res_filt <- all_interval_data %>% select("seqid",
                                            "average_r2_ld_282",
                                            "ge_noperm_mean_gerp",
                                            "max_protein_expression_23tissues",
                                            "max_rna_expression_23tissues",
                                            "atac_seq_peak_count",
                                            "snp_count_goodman_all",
                                            "range_type",
                                            "interval_size",
                                            "residual_goodman_filtered_expression_uniqueCount")

ge_res_filt$range_type <- as.factor(ge_res_filt$range_type)

ge_res_train_filt <- ge_res_filt %>% filter(seqid != 10)
ge_res_test_filt <- ge_res_filt %>% filter(seqid == 10)

# prepare data
gef_res_train_data_mat <- ge_res_train_filt[,-c("seqid",
                                                "residual_goodman_filtered_expression_uniqueCount")] %>% data.matrix()
gef_res_train_label_mat <- ge_res_train_filt[,"residual_goodman_filtered_expression_uniqueCount"] %>% data.matrix()
gef_res_test_data_mat <- ge_res_test_filt[,-c("seqid",
                                              "residual_goodman_filtered_expression_uniqueCount")] %>% data.matrix()
gef_res_test_label_mat <- ge_res_test_filt[,"residual_goodman_filtered_expression_uniqueCount"] %>% data.matrix()
gef_res_dtrain <- xgb.DMatrix(data = gef_res_train_data_mat, label=gef_res_train_label_mat)
gef_res_dtest <- xgb.DMatrix(data = gef_res_test_data_mat, label=gef_res_test_label_mat)
gef_res_watchlist <- list(train=gef_res_dtrain, test=gef_res_dtest)

# run boosting
gef_res_boost <- xgb.train(data=gef_res_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=gef_res_watchlist)
gef_res_pred <- predict(gef_res_boost, gef_res_dtest)
gef_res_importance <- xgb.importance(model = gef_res_boost)


# clean up variables
rm(ge_res_filt)
rm(gef_res_train_data_mat)
rm(gef_res_train_label_mat)
rm(gef_res_test_data_mat)
rm(ge_res_train_filt)
rm(ge_res_test_filt)
rm(gef_res_boost)
rm(gef_res_watchlist)
rm(gef_res_dtrain)
rm(gef_res_dtest)

# ------------------------------
# Goodman Field adjusted (residual) counts, permuted
# ------------------------------
ge_res_permuted <- goodman_expression_res %>% select("seqid",
                                                     "range_type", 
                                                     "interval_size",
                                                     "average_r2_ld_282",
                                                     "snp_count_goodman_all",
                                                     "ge_noperm_mean_gerp",
                                                     "max_protein_expression_23tissues",
                                                     "max_rna_expression_23tissues",
                                                     "atac_seq_peak_count",
                                                     "pleiotropy_gpe")


ge_res_permuted$range_type <- as.factor(ge_res_permuted$range_type)

ge_res_train_permuted <- ge_res_permuted %>% filter(seqid != 10)
ge_res_test_permuted <- ge_res_permuted %>% filter(seqid == 10)

# prepare data
ger_res_train_data_mat <- ge_res_train_permuted %>%
  select(-c("seqid", "pleiotropy_gpe")) %>% data.matrix()
ger_res_train_label_mat <- ge_res_train_permuted[,"pleiotropy_gpe"] %>% data.matrix()
ger_res_test_data_mat <- ge_res_test_permuted %>%
  select(-c("seqid", "pleiotropy_gpe")) %>% data.matrix()
ger_res_test_label_mat <- ge_res_test_permuted[,"pleiotropy_gpe"] %>% data.matrix()
ger_res_dtrain <- xgb.DMatrix(data = ger_res_train_data_mat, label=ger_res_train_label_mat)
ger_res_dtest <- xgb.DMatrix(data = ger_res_test_data_mat, label=ger_res_test_label_mat)
ger_res_watchlist <- list(train=ger_res_dtrain, test=ger_res_dtest)

# run boosting
ger_res_boost <- xgb.train(data=ger_res_dtrain, max.depth=6, eta=0.5, nthread = 2, nrounds=8, watchlist=ger_res_watchlist)
ger_res_pred <- predict(ger_res_boost, ger_res_dtest)
ger_res_importance <- xgb.importance(model = ger_res_boost)


# clean up variables
rm(ge_res_permuted)
rm(ger_res_train_data_mat)
rm(ger_res_train_label_mat)
rm(ger_res_test_data_mat)
rm(ge_res_train_permuted)
rm(ge_res_test_permuted)
rm(ger_res_boost)
rm(ger_res_watchlist)
rm(ger_res_dtrain)
rm(ger_res_dtest)


# --------------------------------------------------------------------------------
# Make a table summarizing results
# --------------------------------------------------------------------------------

# Returns MSE and R^2
fit.stats <- function(pred, labels) {
  c(mean((pred - labels)^2), cor(pred, labels)^2) %>% return
}

adj <- data.frame(matrix(ncol = 3, nrow = 8))
colnames(adj) <- c("model", "MSE", "r_squared")
adj$model <- models
adj[1,2:3] <- fit.stats(npf_res_pred, npf_res_test_label_mat)
adj[2,2:3] <- fit.stats(npr_res_pred, npr_res_test_label_mat)
adj[3,2:3] <- fit.stats(gpf_res_pred, gpf_res_test_label_mat)
adj[4,2:3] <- fit.stats(gpr_res_pred, gpr_res_test_label_mat)
adj[5,2:3] <- fit.stats(gmf_res_pred, gmf_res_test_label_mat)
adj[6,2:3] <- fit.stats(gmr_res_pred, gmr_res_test_label_mat)
adj[7,2:3] <- fit.stats(gef_res_pred, gef_res_test_label_mat)
adj[8,2:3] <- fit.stats(ger_res_pred, ger_res_test_label_mat)
data.table::fwrite(adj, "xgboost_prediction_accuracy_adjusted.csv")

