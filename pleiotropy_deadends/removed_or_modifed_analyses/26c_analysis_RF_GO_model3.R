# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-03-13
# Updated... 2022-05-30
#
# Description: Model 3 (RF with GO)
# Using pleiotropy scores from number of unique traits
# combined with many responses in a random forest model
# 
# Running 8 models, 4: Observed data, 4: permuted data
# Each of the 4 models comes from the 4 population-trait categories
#
# Model:
# pleiotropy score ~ ATACseq count + max RNA + max protein + avg. GERP
#                   + interval type + average LD + interval size 
#                   + input SNP count + adjusted expression pleiotropy
#                   + GO terms for specific hypotheses
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


# -------------------------------------------------------------------------------
#                   Add gene ontology terms to RF models
# -------------------------------------------------------------------------------

# Load in subset of go terms that I want to test
sub_go <- read.csv("~/git_projects/pleiotropy/data/go_terms_rf.csv")
sub_go$desc <- paste0(sub_go$desc, " (", sub_go$ontology, ")") # add BP or MF tag to names

# Load in maize gene ontology terms 
# system("wget https://download.maizegdb.org/GeneFunction_and_Expression/Pannzer_GO_Terms/B73%20_v4_GO.out ~/Downloads/")
go <- data.table::fread("~/Downloads/B73 _v4_GO.out")
go$v4_gene <- gsub("_.*", "", go$qpid) # Parse out gene names

# Parse out protein/transcript number
go$trans <- gsub(".*_P", "", go$qpid) %>% as.character() %>% as.numeric()
go <- go %>% filter(trans == 1) %>% select(v4_gene, goid, ontology, desc)
go$goid <- paste0("GO:", stringr::str_pad(go$goid, 7, pad = "0")) # Turn go numbers back into GO IDS

# Merging together the two files to get more metadata (i.e. what genes the selected GO terms match up with)
parent_go <- merge(go, sub_go, by = "goid")
parent_go <- parent_go[,c(1,2)] %>% unique()

# Create a 'seen' column 
# (1 = that genic interval matches that term, 0 = that genic and all intergenic intervals do not have that term)
parent_go$seen <- rep(1, nrow(parent_go))
parent_go$goid <- gsub(":", "_", parent_go$goid)

# cast data and make wider
parent_go <- tidyr::pivot_wider(parent_go, names_from = goid, values_from = seen, values_fill = 0) 

# read in data aggregatied by aggregate_interval.data.R
# data_dir <- "/Volumes/merrittData1/pleiotropy/interval_data"
data_dir <- "~/git_projects/pleiotropy/data/"
all_interval_data <- data.table::fread(paste0(data_dir, "/all_interval_data_unique_filter.csv"))

# Replace NAs with -1 for GERP columns
all_interval_data$np_noperm_mean_gerp[is.na(all_interval_data$np_noperm_mean_gerp)] <- -1
all_interval_data$gp_noperm_mean_gerp[is.na(all_interval_data$gp_noperm_mean_gerp)] <- -1
all_interval_data$gm_noperm_mean_gerp[is.na(all_interval_data$gm_noperm_mean_gerp)] <- -1
all_interval_data$ge_noperm_mean_gerp[is.na(all_interval_data$ge_noperm_mean_gerp)] <- -1

# Merge ontology terms with interval data, 7k terms
all_interval_data <- merge(all_interval_data, parent_go, by.x = "v4_gene_rna", by.y = "v4_gene", all.x = TRUE)

# If NA in the columns with GO terms, add zero
all_interval_data[, 190:203][is.na(all_interval_data[, 190:203])] <- 0


# -------------------------------------------------------------------
# Plot correlation of predictors
# -------------------------------------------------------------------

# complete correlation matrix ---------------------------------------

# Subset df to just terms in RF models
temp <- all_interval_data %>% select("average_r2_ld_nam",
                                     "average_r2_ld_282",
                                     "snp_count_goodman_all",
                                     "snp_count_nam_all",
                                     "interval_size",
                                     "np_noperm_mean_gerp",
                                     "gp_noperm_mean_gerp",
                                     "gm_noperm_mean_gerp",
                                     "ge_noperm_mean_gerp",
                                     "max_protein_expression_23tissues",
                                     "max_rna_expression_23tissues",
                                     "atac_seq_peak_count",
                                     "nam_filtered_all_true_uniqueCount", "npp_unique_perm_1",
                                     "goodman_filtered_physiological_true_uniqueCount", "gpp_unique_perm_1",
                                     "goodman_filtered_metabolite_true_uniqueCount","gpm_perm_1",
                                     "goodman_filtered_expression_uniqueCount","gpe_perm_1",
                                     "residual_goodman_filtered_expression_uniqueCount",
                                     "GO_0000166","GO_0003677","GO_0003700","GO_0003723","GO_0006351",
                                     "GO_0006355","GO_0006412","GO_0006417","GO_0007165","GO_0009908",
                                     "GO_0009909","GO_0048366","GO_0050790","GO_2000024")
# Create correlations
corr <- round(cor(temp), 1)

# Change rownames to look nice in plot
rename_corr <- c("Average LD NAM",
                 "Average LD GAP",
                 "Input SNPs GAP",
                 "Input SNPs NAM",
                 "Interval Size",
                 "GERP NAM Field",
                 "GERP GAP Field",
                 "GERP GAP Mass Features",
                 "GERP GAP Expression",
                 "Max Protein Expression",
                 "Max RNA Expression",
                 "ATACseq Peak Count",
                 "NAM Field Pleiotropy",
                 "NAM Field Perm 1 Pleiotropy",
                 "GAP Field Pleiotropy",
                 "GAP Field Perm 1 Pleiotropy",
                 "GAP Mass Feature Pleiotropy",
                 "GAP Mass Feature Perm 1 Pleiotropy",
                 "GAP Expression Pleiotropy",
                 "GAP Expression Perm 1 Pleiotropy",
                 "Adj. GAP Expression Pleiotropy",
                 "GO_0000166","GO_0003677","GO_0003700","GO_0003723","GO_0006351",
                 "GO_0006355","GO_0006412","GO_0006417","GO_0007165","GO_0009908",
                 "GO_0009909","GO_0048366","GO_0050790","GO_2000024")
colnames(corr) <- rename_corr
rownames(corr) <- rename_corr

# plot plot
a <- ggcorrplot::ggcorrplot(corr, hc.order = FALSE, type = "lower",
                            lab = TRUE,
                            outline.col = "black",
                            ggtheme = ggplot2::theme_minimal,
                            colors = c("#6D9EC1", "white", "#E46726"))

ggsave("~/git_projects/pleiotropy/images/correlation_matrix_rf_terms.png",
       plot = a,
       width = 17,
       height = 13, 
       units = "in",
       dpi = "retina")

# subsampled correlation matrix

# Subset df to just terms in RF models (minus GO terms)
temp <- all_interval_data %>% select("average_r2_ld_nam",
                                     "average_r2_ld_282",
                                     "snp_count_goodman_all",
                                     "snp_count_nam_all",
                                     "interval_size",
                                     "np_noperm_mean_gerp",
                                     "gp_noperm_mean_gerp",
                                     "gm_noperm_mean_gerp",
                                     "ge_noperm_mean_gerp",
                                     "max_protein_expression_23tissues",
                                     "max_rna_expression_23tissues",
                                     "atac_seq_peak_count",
                                     "nam_filtered_all_true_uniqueCount", "npp_unique_perm_1",
                                     "goodman_filtered_physiological_true_uniqueCount", "gpp_unique_perm_1",
                                     "goodman_filtered_metabolite_true_uniqueCount","gpm_perm_1",
                                     "goodman_filtered_expression_uniqueCount","gpe_perm_1",
                                     "residual_goodman_filtered_expression_uniqueCount")
# Create correlations
corr <- round(cor(temp), 1)

# Change rownames to look nice in plot
rename_corr <- c("Average LD NAM",
                 "Average LD GAP",
                 "Input SNPs GAP",
                 "Input SNPs NAM",
                 "Interval Size",
                 "GERP NAM Field",
                 "GERP GAP Field",
                 "GERP GAP Mass Features",
                 "GERP GAP Expression",
                 "Max Protein Expression",
                 "Max RNA Expression",
                 "ATACseq Peak Count",
                 "NAM Field Pleiotropy",
                 "NAM Field Perm 1 Pleiotropy",
                 "GAP Field Pleiotropy",
                 "GAP Field Perm 1 Pleiotropy",
                 "GAP Mass Features Pleiotropy",
                 "GAP Mass Features Perm 1 Pleiotropy",
                 "GAP Expression Pleiotropy",
                 "GAP Expression Perm 1 Pleiotropy",
                 "Adj. GAP Expression Pleiotropy")
colnames(corr) <- rename_corr
rownames(corr) <- rename_corr

# plot plot
b <- ggcorrplot::ggcorrplot(corr, hc.order = FALSE, type = "lower",
                            lab = TRUE,
                            outline.col = "black",
                            ggtheme = ggplot2::theme_minimal,
                            colors = c("#6D9EC1", "white", "#E46726"))

ggsave("~/git_projects/pleiotropy/images/correlation_matrix_rf_terms_subsampled.png",
       plot = b,
       width = 13,
       height = 9, 
       units = "in",
       dpi = "retina")


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
         "np_noperm_mean_gerp",
         "gp_noperm_mean_gerp",
         "gm_noperm_mean_gerp",
         "ge_noperm_mean_gerp",
         "max_protein_expression_23tissues",
         "max_rna_expression_23tissues",
         "atac_seq_peak_count", "residual_goodman_filtered_expression_uniqueCount",
         "GO_0000166","GO_0003677","GO_0003700","GO_0003723","GO_0006351","GO_0006355",
         "GO_0006412","GO_0006417","GO_0007165","GO_0009908","GO_0009909","GO_0048366",
         "GO_0050790","GO_2000024",
         "npp_unique_perm_1","npp_unique_perm_2","npp_unique_perm_3","npp_unique_perm_4","npp_unique_perm_5",
         "npp_unique_perm_6","npp_unique_perm_7","npp_unique_perm_8","npp_unique_perm_9","npp_unique_perm_10",
         "gpp_unique_perm_1","gpp_unique_perm_2","gpp_unique_perm_3","gpp_unique_perm_4","gpp_unique_perm_5",
         "gpp_unique_perm_6","gpp_unique_perm_7","gpp_unique_perm_8","gpp_unique_perm_9","gpp_unique_perm_10",
         "gpm_perm_1","gpm_perm_2","gpm_perm_3","gpm_perm_4","gpm_perm_5",
         "gpm_perm_6","gpm_perm_7","gpm_perm_8","gpm_perm_9","gpm_perm_10",
         "gpe_perm_1","gpe_perm_2","gpe_perm_3","gpe_perm_4","gpe_perm_5")

# Melt individual datasets, then merge 
# (because I can't figure out how to wrap this all in a single pivot_longer command)
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
                          list(nam_field[,c(1:30,58,57)], goodman_field[,c(1:2,57:58)], goodman_mass[,c(1:2,57:58)]))
dim(permuted_merged) # should be: 754900     34
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

np_obs <- all_interval_data %>% select("rr_id", "seqid", 
                                       "average_r2_ld_nam", 
                                       "np_noperm_mean_gerp",
                                       "max_protein_expression_23tissues",
                                       "max_rna_expression_23tissues", 
                                       "atac_seq_peak_count",
                                       "snp_count_nam_all", 
                                       "range_type", "interval_size",
                                       "residual_goodman_filtered_expression_uniqueCount",
                                       "nam_filtered_all_true_uniqueCount",
                                       "GO_0000166","GO_0003677","GO_0003700","GO_0003723",
                                       "GO_0006351","GO_0006355","GO_0006412","GO_0006417",
                                       "GO_0007165","GO_0009908","GO_0009909","GO_0048366",
                                       "GO_0050790","GO_2000024")
colnames(np_obs)
dim(np_obs) # 75490    26
np_obs$range_type <- as.factor(np_obs$range_type)

# ------------------------------
# NAM Permuted
# ------------------------------

np_permuted <- permuted_merged %>% select("rr_id", "range_type" , "seqid",
                                          "interval_size", "average_r2_ld_nam",
                                          "snp_count_nam_all", "np_noperm_mean_gerp",
                                          "max_protein_expression_23tissues",
                                          "max_rna_expression_23tissues",
                                          "atac_seq_peak_count",
                                          "residual_goodman_filtered_expression_uniqueCount",
                                          "pleiotropy_npp",
                                          "GO_0000166","GO_0003677","GO_0003700","GO_0003723",
                                          "GO_0006351","GO_0006355","GO_0006412","GO_0006417",
                                          "GO_0007165","GO_0009908","GO_0009909","GO_0048366",
                                          "GO_0050790","GO_2000024")
colnames(np_permuted)
dim(np_permuted) # 754900     26
np_permuted$range_type <- as.factor(np_permuted$range_type)


# ------------------------------
# goodman Field Observed Data
# ------------------------------

gp_obs <- all_interval_data %>% select("rr_id", "seqid", 
                                       "average_r2_ld_282", 
                                       "gp_noperm_mean_gerp",
                                       "max_protein_expression_23tissues",
                                       "max_rna_expression_23tissues", 
                                       "atac_seq_peak_count",
                                       "snp_count_goodman_all", 
                                       "range_type", "interval_size",
                                       "residual_goodman_filtered_expression_uniqueCount",
                                       "goodman_filtered_physiological_true_uniqueCount",
                                       "GO_0000166","GO_0003677","GO_0003700","GO_0003723",
                                       "GO_0006351","GO_0006355","GO_0006412","GO_0006417",
                                       "GO_0007165","GO_0009908","GO_0009909","GO_0048366",
                                       "GO_0050790","GO_2000024")
colnames(gp_obs)
dim(gp_obs) # 75490    26
gp_obs$range_type <- as.factor(gp_obs$range_type)

# ------------------------------
# Goodman Field Permuted
# ------------------------------

gp_permuted <- permuted_merged %>% select("rr_id", "range_type" , "seqid",
                                          "interval_size", "average_r2_ld_282",
                                          "snp_count_goodman_all", "gp_noperm_mean_gerp",
                                          "max_protein_expression_23tissues",
                                          "max_rna_expression_23tissues",
                                          "atac_seq_peak_count",
                                          "residual_goodman_filtered_expression_uniqueCount",
                                          "pleiotropy_gpp",
                                          "GO_0000166","GO_0003677","GO_0003700","GO_0003723",
                                          "GO_0006351","GO_0006355","GO_0006412","GO_0006417",
                                          "GO_0007165","GO_0009908","GO_0009909","GO_0048366",
                                          "GO_0050790","GO_2000024")
colnames(gp_permuted)
dim(gp_permuted) # 754900     26
gp_permuted$range_type <- as.factor(gp_permuted$range_type)


# ------------------------------
# goodman Mass Features Observed Data
# ------------------------------

gm_obs <- all_interval_data %>% select("rr_id", "seqid", 
                                       "average_r2_ld_282", 
                                       "gm_noperm_mean_gerp",
                                       "max_protein_expression_23tissues",
                                       "max_rna_expression_23tissues", 
                                       "atac_seq_peak_count",
                                       "snp_count_goodman_all", 
                                       "range_type", "interval_size",
                                       "residual_goodman_filtered_expression_uniqueCount",
                                       "goodman_filtered_metabolite_true_uniqueCount",
                                       "GO_0000166","GO_0003677","GO_0003700","GO_0003723",
                                       "GO_0006351","GO_0006355","GO_0006412","GO_0006417",
                                       "GO_0007165","GO_0009908","GO_0009909","GO_0048366",
                                       "GO_0050790","GO_2000024")
colnames(gm_obs)
dim(gm_obs)# 75490    26
gm_obs$range_type <- as.factor(gm_obs$range_type)


# ------------------------------
# Goodman Mass Features Permuted
# ------------------------------

gm_permuted <- permuted_merged %>% select("rr_id", "range_type" , "seqid",
                                          "interval_size", "average_r2_ld_282",
                                          "snp_count_goodman_all", "gm_noperm_mean_gerp",
                                          "max_protein_expression_23tissues",
                                          "max_rna_expression_23tissues",
                                          "atac_seq_peak_count",
                                          "residual_goodman_filtered_expression_uniqueCount",
                                          "pleiotropy_gpm",
                                          "GO_0000166","GO_0003677","GO_0003700","GO_0003723",
                                          "GO_0006351","GO_0006355","GO_0006412","GO_0006417",
                                          "GO_0007165","GO_0009908","GO_0009909","GO_0048366",
                                          "GO_0050790","GO_2000024")
colnames(gm_permuted)
dim(gm_permuted) # 754900     26
gm_permuted$range_type <- as.factor(gm_permuted$range_type)


# ------------------------------
# goodman expression Observed Data
# ------------------------------

ge_obs <- all_interval_data %>% select("rr_id", "seqid", 
                                       "average_r2_ld_282", 
                                       "ge_noperm_mean_gerp",
                                       "max_protein_expression_23tissues",
                                       "max_rna_expression_23tissues", 
                                       "atac_seq_peak_count",
                                       "snp_count_goodman_all", 
                                       "range_type", "interval_size",
                                       "goodman_filtered_expression_uniqueCount",
                                       "GO_0000166","GO_0003677","GO_0003700","GO_0003723",
                                       "GO_0006351","GO_0006355","GO_0006412","GO_0006417",
                                       "GO_0007165","GO_0009908","GO_0009909","GO_0048366",
                                       "GO_0050790","GO_2000024")
colnames(ge_obs)
dim(ge_obs) # 75490    25
ge_obs$range_type <- as.factor(ge_obs$range_type)


# ------------------------------
# Goodman Expression Permuted
# ------------------------------

ge_permuted <- goodman_expression %>% select("rr_id", 
                                             "range_type", 
                                             "seqid", 
                                             "average_r2_ld_282",
                                             "interval_size",
                                             "snp_count_goodman_all",
                                             "gp_noperm_mean_gerp", 
                                             "max_protein_expression_23tissues",
                                             "max_rna_expression_23tissues",
                                             "atac_seq_peak_count",
                                             "pleiotropy_gpe",
                                             "GO_0000166","GO_0003677","GO_0003700","GO_0003723",
                                             "GO_0006351","GO_0006355","GO_0006412","GO_0006417",
                                             "GO_0007165","GO_0009908","GO_0009909","GO_0048366",
                                             "GO_0050790","GO_2000024")
colnames(ge_permuted)
dim(ge_permuted) # 377450     25
ge_permuted$range_type <- as.factor(ge_permuted$range_type)


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

# Load in subset of go terms that I want to test
go_desc <- read.csv("~/git_projects/pleiotropy/data/go_terms_rf.csv")
go_desc$desc <- paste0(go_desc$desc, " (", go_desc$ontology, ")") # add BP or MF tag to names
go_desc$goid <- gsub(":", "_", go_desc$goid)
go_desc <- go_desc %>% select(-"ontology")
colnames(go_desc) <- c("old", "Variable")
go_desc$Type <- rep("Biological", nrow(go_desc))

# Add go term descriptions
pretty_names <- rbind(pretty_names, go_desc)


# ---------------------------------------------------------------------------------
#                       Run Random Forest Models
# ---------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
#                               nam Field 
# --------------------------------------------------------------------------------

# Observed Data model ----------------------------------
nam_field_obs_model3 <- ranger::ranger(nam_filtered_all_true_uniqueCount ~ ., 
                                       importance = 'impurity',
                                       case.weights = np_case_obs,
                                       holdout = TRUE,
                                       data = np_obs[,-c(1,2)],
                                       num.trees = 500)
nam_field_obs_model3

# prep for importance plot
nam_field_obs_model3_importance <- as.data.frame(nam_field_obs_model3$variable.importance)
nam_field_obs_model3_importance$type <- rownames(nam_field_obs_model3_importance)
colnames(nam_field_obs_model3_importance) <- c("Importance", "Variable_Type")
nam_field_obs_model3_importance$Importance <- nam_field_obs_model3_importance$Importance/max(nam_field_obs_model3_importance$Importance)
nam_field_obs_model3_importance <- merge(x = nam_field_obs_model3_importance, y = pretty_names, 
                                         by.x = "Variable_Type", by.y = "old", all.x = TRUE)

# Permuted model-------------------------------------
nam_field_perm_model3 <- ranger::ranger(pleiotropy_npp ~ ., 
                                        importance = 'impurity',
                                        case.weights = np_case_permuted,
                                        holdout = TRUE,
                                        data = np_permuted[,-c(1,3)],
                                        num.trees = 500)
nam_field_perm_model3

# prep for importance plot
nam_field_perm_model3_importance <- as.data.frame(nam_field_perm_model3$variable.importance)
nam_field_perm_model3_importance$type <- rownames(nam_field_perm_model3_importance)
colnames(nam_field_perm_model3_importance) <- c("Importance", "Variable_Type")
nam_field_perm_model3_importance$Importance <- nam_field_perm_model3_importance$Importance/max(nam_field_perm_model3_importance$Importance)
nam_field_perm_model3_importance <- merge(x = nam_field_perm_model3_importance, y = pretty_names, 
                                          by.x = "Variable_Type", by.y = "old", all.x = TRUE)


# --------------------------------------------------------------------------------
#                               282 Field
# --------------------------------------------------------------------------------

# 282 Field Observed Data-----------------------------
goodman_field_obs_model3 <- ranger::ranger(goodman_filtered_physiological_true_uniqueCount ~ ., 
                                           importance = 'impurity',
                                           case.weights = gp_case_obs,
                                           holdout = TRUE,
                                           data = gp_obs[,-c(1,2)],
                                           num.trees = 500)
goodman_field_obs_model3

# prep for importance plot
goodman_field_obs_model3_importance <- as.data.frame(goodman_field_obs_model3$variable.importance)
goodman_field_obs_model3_importance$type <- rownames(goodman_field_obs_model3_importance)
colnames(goodman_field_obs_model3_importance) <- c("Importance", "Variable_Type")
goodman_field_obs_model3_importance$Importance <- goodman_field_obs_model3_importance$Importance/max(goodman_field_obs_model3_importance$Importance)
goodman_field_obs_model3_importance <- merge(x = goodman_field_obs_model3_importance, y = pretty_names, 
                                             by.x = "Variable_Type", by.y = "old", all.x = TRUE)


# 282 Field Permuted----------------------------------
goodman_field_perm_model3 <- ranger::ranger(pleiotropy_gpp ~ ., 
                                            importance = 'impurity',
                                            case.weights = gp_case_permuted,
                                            holdout = TRUE,
                                            data = gp_permuted[,-c(1,3)],
                                            num.trees = 500)
goodman_field_perm_model3

# prep for importance plot
goodman_field_perm_model3_importance <- as.data.frame(goodman_field_perm_model3$variable.importance)
goodman_field_perm_model3_importance$type <- rownames(goodman_field_perm_model3_importance)
colnames(goodman_field_perm_model3_importance) <- c("Importance", "Variable_Type")
goodman_field_perm_model3_importance$Importance <- goodman_field_perm_model3_importance$Importance/max(goodman_field_perm_model3_importance$Importance)
goodman_field_perm_model3_importance <- merge(x = goodman_field_perm_model3_importance, y = pretty_names, 
                                              by.x = "Variable_Type", by.y = "old", all.x = TRUE)


# --------------------------------------------------------------------------------
#                               282 Mass Features
# --------------------------------------------------------------------------------

# Observed Data model-----------------------------------
goodman_mass_obs_model3 <- ranger::ranger(goodman_filtered_metabolite_true_uniqueCount ~ ., 
                                          importance = 'impurity',
                                          case.weights = gm_case_obs,
                                          holdout = TRUE,
                                          data = gm_obs[,-c(1,2)],
                                          num.trees = 500)
goodman_mass_obs_model3

# prep for importance plot
goodman_mass_obs_model3_importance <- as.data.frame(goodman_mass_obs_model3$variable.importance)
goodman_mass_obs_model3_importance$type <- rownames(goodman_mass_obs_model3_importance)
colnames(goodman_mass_obs_model3_importance) <- c("Importance", "Variable_Type")
goodman_mass_obs_model3_importance$Importance <- goodman_mass_obs_model3_importance$Importance/max(goodman_mass_obs_model3_importance$Importance)
goodman_mass_obs_model3_importance <- merge(x = goodman_mass_obs_model3_importance, y = pretty_names, 
                                            by.x = "Variable_Type", by.y = "old", all.x = TRUE)

# Permuted model-----------------------------------------
goodman_mass_perm_model3 <- ranger::ranger(pleiotropy_gpm ~ ., 
                                           importance = 'impurity',
                                           case.weights = gm_case_permuted,
                                           holdout = TRUE,
                                           data = gm_permuted[,-c(1,3)],
                                           num.trees = 500)
goodman_mass_perm_model3

# prep for importance plot
goodman_mass_perm_model3_importance <- as.data.frame(goodman_mass_perm_model3$variable.importance)
goodman_mass_perm_model3_importance$type <- rownames(goodman_mass_perm_model3_importance)
colnames(goodman_mass_perm_model3_importance) <- c("Importance", "Variable_Type")
goodman_mass_perm_model3_importance$Importance <- goodman_mass_perm_model3_importance$Importance/max(goodman_mass_perm_model3_importance$Importance)
goodman_mass_perm_model3_importance <- merge(x = goodman_mass_perm_model3_importance, y = pretty_names, 
                                             by.x = "Variable_Type", by.y = "old", all.x = TRUE)


# --------------------------------------------------------------------------------
#                                 282 expression
# --------------------------------------------------------------------------------

# Observed Data model
goodman_expr_obs_model3 <- ranger::ranger(goodman_filtered_expression_uniqueCount ~ ., 
                                          importance = 'impurity',
                                          case.weights = ge_case_obs,
                                          holdout = TRUE,
                                          data = ge_obs[,-c(1,2)],
                                          num.trees = 500)
goodman_expr_obs_model3


# prep for importance plot
goodman_expr_obs_model3_importance <- as.data.frame(goodman_expr_obs_model3$variable.importance)
goodman_expr_obs_model3_importance$type <- rownames(goodman_expr_obs_model3_importance)
colnames(goodman_expr_obs_model3_importance) <- c("Importance", "Variable_Type")
goodman_expr_obs_model3_importance$Importance <- goodman_expr_obs_model3_importance$Importance/max(goodman_expr_obs_model3_importance$Importance)
goodman_expr_obs_model3_importance <- merge(x = goodman_expr_obs_model3_importance, y = pretty_names, 
                                            by.x = "Variable_Type", by.y = "old", all.x = TRUE)

# Permuted model
goodman_expr_perm_model3 <- ranger::ranger(pleiotropy_gpe ~ ., 
                                           importance = 'impurity',
                                           case.weights = ge_case_permuted,
                                           holdout = TRUE,
                                           data = ge_permuted[,-c(1,3)],
                                           num.trees = 500)
goodman_expr_perm_model3


# prep for importance plot
goodman_expr_perm_model3_importance <- as.data.frame(goodman_expr_perm_model3$variable.importance)
goodman_expr_perm_model3_importance$type <- rownames(goodman_expr_perm_model3_importance)
colnames(goodman_expr_perm_model3_importance) <- c("Importance", "Variable_Type")
goodman_expr_perm_model3_importance$Importance <- goodman_expr_perm_model3_importance$Importance/max(goodman_expr_perm_model3_importance$Importance)
goodman_expr_perm_model3_importance <- merge(x = goodman_expr_perm_model3_importance, y = pretty_names, 
                                             by.x = "Variable_Type", by.y = "old", all.x = TRUE)


# --------------------------------------------------------------------------------
# Make a table summarizing results
# --------------------------------------------------------------------------------

nam_field_obs_model3
nam_field_perm_model3
goodman_field_obs_model3
goodman_field_perm_model3
goodman_mass_obs_model3
goodman_mass_perm_model3
goodman_expr_obs_model3
goodman_expr_perm_model3

# Assemble model results
temp <- rbind(nam_field_obs_model3[c("prediction.error", "r.squared")] %>% data.frame(), 
              nam_field_perm_model3[c("prediction.error", "r.squared")] %>% data.frame(),
              goodman_field_obs_model3[c("prediction.error", "r.squared")] %>% data.frame(),
              goodman_field_perm_model3[c("prediction.error", "r.squared")] %>% data.frame(),
              goodman_mass_obs_model3[c("prediction.error", "r.squared")] %>% data.frame(),
              goodman_mass_perm_model3[c("prediction.error", "r.squared")] %>% data.frame(),
              goodman_expr_obs_model3[c("prediction.error", "r.squared")] %>% data.frame(),
              goodman_expr_perm_model3[c("prediction.error", "r.squared")] %>% data.frame())
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
data.table::fwrite(temp, "~/git_projects/pleiotropy/data/model3_RF_GO_prediction_accuracy_.csv")


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
a_model3 <- nam_field_obs_model3_importance %>% 
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
b_model3 <- nam_field_perm_model3_importance %>% 
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
  labs(tag = "B")

# goodman Field Observed Data ----------------------------------
c_model3 <- goodman_field_obs_model3_importance %>% 
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
  labs(tag = "C")

# goodman Field Permuted ------------------------------------
d_model3 <- goodman_field_perm_model3_importance %>% 
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
  labs(tag = "D")

# goodman Mass Features Observed Data -------------------------------------
e_model3 <- goodman_mass_obs_model3_importance %>% 
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
  labs(tag = "E")

# goodman Mass Features Permuted --------------------------------------
f_model3 <- goodman_mass_perm_model3_importance %>% 
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
  labs(tag = "F")

# goodman expression Observed Data ------------------------------------
g_model3 <- goodman_expr_obs_model3_importance %>% 
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
  labs(tag = "G")

# goodman expression Permuted --------------------------------------
h_model3 <- goodman_expr_perm_model3_importance %>% 
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
  labs(tag = "H")


# test combined plot ----------------------------------------------
# ggpubr::ggarrange(a_model3,b_model3,c_model3,d_model3,
#                   e_model3,f_model3,g_model3,h_model3,
#                   nrow = 4, ncol = 2,
#                   common.legend = TRUE, legend = "bottom", align = "v")

# Export to file---------------------------------------------------
ggsave("~/git_projects/pleiotropy/images/model3_random_forest_importance_observed_vs_permuted.png",
       plot = ggpubr::ggarrange(a_model3,b_model3,c_model3,d_model3,
                                e_model3,f_model3,g_model3,h_model3,
                                nrow = 4, ncol = 2, common.legend = TRUE, legend = "bottom", align = "v"),
       width = 18,
       height = 25,
       units = "in",
       dpi = "retina")

# Observed data
ggsave("~/git_projects/pleiotropy/images/model3_random_forest_importance_observed.png",
       plot = ggpubr::ggarrange(a_model3,c_model3,e_model3,g_model3,
                                nrow = 2, ncol = 2,
                                common.legend = TRUE, legend = "bottom", align = "v"),
       width = 18,
       height = 14,
       units = "in",
       dpi = "retina")

# Permuted data
ggsave("~/git_projects/pleiotropy/images/model3_random_forest_importance_permuted.png",
       plot = ggpubr::ggarrange(b_model3,d_model3,f_model3,h_model3,
                                nrow = 2, ncol = 2,
                                common.legend = TRUE, legend = "bottom", align = "v"),
       width = 18,
       height = 14,
       units = "in",
       dpi = "retina")


# --------------------------------------------------------------------------------
# Plot & combine all model performance plots 
# --------------------------------------------------------------------------------

# nam phys Observed Data --------------------------------------------------------------
np_pred_obs <- predict(object = nam_field_obs_model3, data = np_test_obs)
np_test_obs$np_pred_obs <- np_pred_obs$predictions

test2 <- paste("~R^2==~", round(cor(np_test_obs$nam_filtered_all_true_uniqueCount, np_test_obs$np_pred_obs)^2, digits = 3))

a_pred_model3 <- ggplot(np_test_obs, aes(x = nam_filtered_all_true_uniqueCount, y = np_pred_obs)) +
  geom_point(colour = "#999999") +
  geom_smooth(method = lm, color = "black", linetype = "solid") + 
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  labs(x = "Observed Pleiotropy", 
       y = "Predicted Pleiotropy", 
       title = "NAM Field Observed") +
  annotate(geom="text", label = test2, x = r2_text_x, y = 70, size = annotate_size, parse = TRUE, hjust = 0) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  ylim(0,90) +
  xlim(0,90) +
  labs(tag = "A")


# nam phys Permuted-----------------------------------------------------------------
np_pred_permuted <- predict(object = nam_field_perm_model3, data = np_test_permuted)
np_test_permuted$np_pred_permuted <- np_pred_permuted$predictions

test2 <- paste("~R^2==~", 
               round(cor(np_test_permuted$pleiotropy_npp, np_test_permuted$np_pred_permuted)^2, digits = 3))

b_pred_model3 <- ggplot(np_test_permuted, aes(x = pleiotropy_npp, y = np_pred_permuted)) +
  geom_point(colour = "#999999") +
  geom_smooth(method = lm, color = "black") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Permuted Pleiotropy", 
       y = "Predicted Pleiotropy", 
       title = "NAM Field Permuted") +
  annotate(geom="text", label = test2, x = r2_text_x, y = 70, size = annotate_size, parse = TRUE, hjust = 0) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  ylim(0,90) +
  xlim(0,90) +
  labs(tag = "B")

# goodman Field Observed Data--------------------------------------------------
gp_pred_obs <- predict(object = goodman_field_obs_model3, data = gp_test_obs)
gp_test_obs$gp_pred_obs <- gp_pred_obs$predictions

test2 <- paste("~R^2==~", 
               round(cor(gp_test_obs$goodman_filtered_physiological_true_uniqueCount, gp_test_obs$gp_pred_obs)^2, digits = 3))

c_pred_model3 <- ggplot(gp_test_obs, aes(x = goodman_filtered_physiological_true_uniqueCount, y = gp_pred_obs)) +
  geom_point(colour = "#999999") +
  geom_smooth(method = lm, color = "black") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Observed Pleiotropy", 
       y = "Predicted Pleiotropy", 
       title = "GAP Field Observed") +
  annotate(geom="text", label = test2, x = r2_text_x, y = 65, size = annotate_size, parse = TRUE, hjust = 0) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  ylim(0,85) +
  xlim(0,85) +
  labs(tag = "C")

# goodman Field Permuted---------------------------------------------------
gp_pred_permuted <- predict(object = goodman_field_perm_model3, data = gp_test_permuted)
gp_test_permuted$gp_pred_permuted <- gp_pred_permuted$predictions

test2 <- paste("~R^2==~", 
               round(cor(gp_test_permuted$pleiotropy_gpp, gp_test_permuted$gp_pred_permuted)^2, digits = 3))

d_pred_model3 <- ggplot(gp_test_permuted, aes(x = pleiotropy_gpp, y = gp_pred_permuted)) +
  geom_point(colour = "#999999") +
  geom_smooth(method = lm, color = "black") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Permuted Pleiotropy", 
       y = "Predicted Pleiotropy", 
       title = "GAP Field Permuted") +
  annotate(geom="text", label = test2, x = r2_text_x, y = 65, size = annotate_size, parse = TRUE, hjust = 0) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  ylim(0,85) +
  xlim(0,85) +
  labs(tag = "D")

# goodman Mass Features Observed Data----------------------------------------------------
# Predictions
gm_pred_obs <- predict(object = goodman_mass_obs_model3, data = gm_test_obs)
gm_test_obs$gm_pred_obs <- gm_pred_obs$predictions

test2 <- paste("~R^2==~", 
               round(cor(gm_test_obs$goodman_filtered_metabolite_true_uniqueCount, gm_test_obs$gm_pred_obs)^2, digits = 3))

e_pred_model3 <- ggplot(gm_test_obs, aes(x = goodman_filtered_metabolite_true_uniqueCount, y = gm_pred_obs)) +
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
  labs(tag = "E")

# goodman Mass Features Permuted-----------------------------------------------------
gm_pred_permuted <- predict(object = goodman_mass_perm_model3, data = gm_test_permuted)
gm_test_permuted$gm_pred_permuted <- gm_pred_permuted$predictions

test2 <- paste("~R^2==~", 
               round(cor(gm_test_permuted$pleiotropy_gpm, gm_test_permuted$gm_pred_permuted)^2, digits = 3))

f_pred_model3 <- ggplot(gm_test_permuted, aes(x = pleiotropy_gpm, y = gm_pred_permuted)) +
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
  labs(tag = "F")

# goodman expression Observed Data---------------------------------------------------
ge_pred_obs <- predict(object = goodman_expr_obs_model3, data = ge_test_obs)
ge_test_obs$ge_pred_obs <- ge_pred_obs$predictions

test2 <- paste("~R^2==~", 
               round(cor(ge_test_obs$goodman_filtered_expression_uniqueCount, ge_test_obs$ge_pred_obs)^2, digits = 3))

g_pred_model3 <- ggplot(ge_test_obs, aes(x = goodman_filtered_expression_uniqueCount, y = ge_pred_obs)) +
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
  labs(tag = "G")

# goodman expression Permuted------------------------------------------------------
ge_pred_permuted <- predict(object = goodman_expr_perm_model3, data = ge_test_permuted)
ge_test_permuted$ge_pred_permuted <- ge_pred_permuted$predictions

test2 <- paste("~R^2==~", 
               round(cor(ge_test_permuted$pleiotropy_gpe, ge_test_permuted$ge_pred_permuted)^2, digits = 3))

h_pred_model3 <- ggplot(ge_test_permuted, aes(x = pleiotropy_gpe, y = ge_pred_permuted)) +
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
  labs(tag = "H")

# test combined plot-------------------------------------------------------------
# ggpubr::ggarrange(a_pred_model3,b_pred_model3,c_pred_model3,d_pred_model3,
#                   e_pred_model3,f_pred_model3,g_pred_model3,h_pred_model3,
#                   nrow = 4, ncol = 2,
#                   common.legend = TRUE, legend = "bottom", align = "v")


# Export to file-----------------------------------------------------------------
ggsave("~/git_projects/pleiotropy/images/model3_random_forest_prediction_accuracy_observed_vs_permuted.png",
       plot = ggpubr::ggarrange(a_pred_model3,b_pred_model3,c_pred_model3,d_pred_model3,
                                e_pred_model3,f_pred_model3,g_pred_model3,h_pred_model3,
                                nrow = 4, ncol = 2, 
                                common.legend = TRUE, legend = "bottom", align = "v"),
       width = 11,
       height = 10, 
       units = "in",
       dpi = "retina")
