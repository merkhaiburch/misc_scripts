# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-12-13 
# Updated... 2022-06-24
#
# Description 
# Seeing if different GO terms have associations with pleiotropy
# Source: https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/
# ---------------------------------------------------------------

# Load in packages
library(magrittr)
library(tidyr)
library(ggplot2)
library(topGO)


# -------------------------------
# Data gathering & formatting
# -------------------------------

# read in data aggregatied by aggregate_interval.data.R
data_dir <- "/Volumes/merrittData1/pleiotropy/interval_data"
all_interval_data <- data.table::fread(paste0(data_dir, "/all_interval_data_unique_filter.csv"))

# Gather v5 go ontologies
go_path = "~/git_projects/pleiotropy/data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.interproscan.tsv"

# Load in v4 --> v5 gene ID table
xref <- read.delim("~/git_projects/pleiotropy/data/B73v4_B73v5_liftoff_genemodel_CDS_xref_shortened.txt", header = F)
colnames(xref) <- c("v4", "v5")
xref <- xref[!grepl("chr",xref$v5),]

# Convert my v4 gene IDs to v5
v5_all_interval_data <- merge(x = all_interval_data, y = xref,
                              by.x = "v4_gene_rna", by.y = "v4")
dim(v5_all_interval_data)
length(unique(v5_all_interval_data$v5))

# Create a list where names are genes and content is all gene ontology terms
all_go = readr::read_tsv(go_path, col_names = FALSE) %>%
  dplyr::select(Tx=X1, GO=X14) %>%
  dplyr::mutate(Gene = sub("_T[0-9]{3}", "", Tx)) %>%
  separate_rows(GO, sep="\\|") %>%
  split(.$Gene) %>%
  purrr::map(~ .x %>% drop_na() %>% dplyr::pull(GO) %>% unique())


# ------------------------------------
#   Set global variables
# ------------------------------------

# Set a threshold for adjusted pleiotropy scores to be highly and 
# lowly pleiotropic for GO analysis
p_val_cutoff <- 0.05

# if wanted to do 25th percentiles
high_cutoff <- 0.75
low_cutoff <- 0.25

# If a 90/10 percentile
# high_cutoff <- 0.9 # high = top 90th percentile
# low_cutoff <- 0.1 # low = bottom 10th percentile

# -------------------------------------------------------------------------------------------------
#       Filter out highly and lowly pleiotropy data for GO Ontology Enrichment
# -------------------------------------------------------------------------------------------------

# nam physiological -----------------------------------------------

# just take values above zero + check distribution
temp <-  v5_all_interval_data %>% 
  dplyr::select("residual_nam_filtered_all_true_uniqueCount", "v5") %>% 
  dplyr::filter(residual_nam_filtered_all_true_uniqueCount > 0)
hist(temp$residual_nam_filtered_all_true_uniqueCount)

# Calculate quantile
temp_high <- quantile(temp$residual_nam_filtered_all_true_uniqueCount, high_cutoff)
temp_low <- quantile(temp$residual_nam_filtered_all_true_uniqueCount, low_cutoff)

# Subset out Observed, highly pleiotropic regions, physiological only
np_observed_highly_names <- temp %>% 
  dplyr::select("residual_nam_filtered_all_true_uniqueCount", "v5") %>% 
  dplyr::filter(residual_nam_filtered_all_true_uniqueCount >= temp_high) %>% 
  dplyr::pull(v5) %>% 
  unique()
np_observed_highly_genes <- factor(as.numeric(names(all_go) %in% np_observed_highly_names))
names(np_observed_highly_genes) = names(all_go)

# Subset out observed, lowly pleiotropic regions, physiological only
np_observed_lowly_names <- temp %>% 
  dplyr::select("residual_nam_filtered_all_true_uniqueCount", "v5") %>% 
  dplyr::filter(residual_nam_filtered_all_true_uniqueCount <= temp_low) %>% 
  dplyr::pull(v5) %>% 
  unique()
np_observed_lowly_genes <- factor(as.numeric(names(all_go) %in% np_observed_lowly_names))
names(np_observed_lowly_genes) = names(all_go)


# goodman physiological -----------------------------------------
temp_gf <-  v5_all_interval_data %>% 
  dplyr::select("residual_goodman_filtered_physiological_true_uniqueCount", "v5") %>% 
  dplyr::filter(residual_goodman_filtered_physiological_true_uniqueCount > 0)
hist(temp_gf$residual_goodman_filtered_physiological_true_uniqueCount)

temp_high <- quantile(temp_gf$residual_goodman_filtered_physiological_true_uniqueCount, high_cutoff)
temp_low <- quantile(temp_gf$residual_goodman_filtered_physiological_true_uniqueCount, low_cutoff)

# Subset out observed, highly pleiotropic regions, physiological only
gp_observed_highly_names <- temp_gf %>% 
  dplyr::select("residual_goodman_filtered_physiological_true_uniqueCount", "v5") %>% 
  dplyr::filter(residual_goodman_filtered_physiological_true_uniqueCount >= temp_high) %>% 
  dplyr::pull(v5) %>% 
  unique()
gp_observed_highly_genes <- factor(as.numeric(names(all_go) %in% gp_observed_highly_names))
names(gp_observed_highly_genes) = names(all_go)

# Subset out observed, lowly pleiotropic regions, physiological only
gp_observed_lowly_names <- temp_gf %>% 
  dplyr::select("residual_goodman_filtered_physiological_true_uniqueCount", "v5") %>% 
  dplyr::filter(residual_goodman_filtered_physiological_true_uniqueCount <= temp_low) %>% 
  dplyr::pull(v5) %>% 
  unique()
gp_observed_lowly_genes <- factor(as.numeric(names(all_go) %in% gp_observed_lowly_names))
names(gp_observed_lowly_genes) = names(all_go)


# Metabolite ----------------------------------------------------
# just take values above zero + check distribution
temp_gm <-  v5_all_interval_data %>% 
  dplyr::select("residual_goodman_filtered_metabolite_true_uniqueCount", "v5") %>% 
  dplyr::filter(residual_goodman_filtered_metabolite_true_uniqueCount > 0)
hist(temp_gf$residual_goodman_filtered_metabolite_true_uniqueCount)

temp_high <- quantile(temp_gm$residual_goodman_filtered_metabolite_true_uniqueCount, high_cutoff)
temp_low <- quantile(temp_gm$residual_goodman_filtered_metabolite_true_uniqueCount, low_cutoff)

# Subset out observed, highly pleiotropic regions, metabolite only
gm_observed_highly_names <- temp_gm %>% 
  dplyr::select("residual_goodman_filtered_metabolite_true_uniqueCount", "v5") %>% 
  dplyr::filter(residual_goodman_filtered_metabolite_true_uniqueCount >= temp_high) %>% 
  dplyr::pull(v5) %>% 
  unique()
gm_observed_highly_genes <- factor(as.numeric(names(all_go) %in% gm_observed_highly_names))
names(gm_observed_highly_genes) = names(all_go)

# Subset out observed, lowly pleiotropic regions, metabolite only
gm_observed_lowly_names <- temp_gm %>% 
  dplyr::select("residual_goodman_filtered_metabolite_true_uniqueCount", "v5") %>% 
  dplyr::filter(residual_goodman_filtered_metabolite_true_uniqueCount <= temp_low) %>% 
  dplyr::pull(v5) %>% 
  unique()
gm_observed_lowly_genes <- factor(as.numeric(names(all_go) %in% gm_observed_lowly_names))
names(gm_observed_lowly_genes) = names(all_go)


# Expression --------------------------------------------------
temp_ge <-  v5_all_interval_data %>% 
  dplyr::select("residual_goodman_filtered_expression_uniqueCount", "v5") %>% 
  dplyr::filter(residual_goodman_filtered_expression_uniqueCount > 0)
hist(temp_ge$residual_goodman_filtered_expression_uniqueCount)

temp_high <- quantile(temp_ge$residual_goodman_filtered_expression_uniqueCount, high_cutoff)
temp_low <- quantile(temp_ge$residual_goodman_filtered_expression_uniqueCount, low_cutoff)

# Subset out observed, highly pleiotropic regions, expression only
ge_observed_highly_names <- temp_ge %>% 
  dplyr::select("residual_goodman_filtered_expression_uniqueCount", "v5") %>% 
  dplyr::filter(residual_goodman_filtered_expression_uniqueCount >= temp_high) %>% 
  dplyr::pull(v5) %>% 
  unique()
ge_observed_highly_genes <- factor(as.numeric(names(all_go) %in% ge_observed_highly_names))
names(ge_observed_highly_genes) = names(all_go)

# Subset out observed, lowly pleiotropic regions, expression only
ge_observed_lowly_names <- temp_ge %>% 
  dplyr::select("residual_goodman_filtered_expression_uniqueCount", "v5") %>% 
  dplyr::filter(residual_goodman_filtered_expression_uniqueCount <= temp_low) %>% 
  dplyr::pull(v5) %>% 
  unique()
ge_observed_lowly_genes <- factor(as.numeric(names(all_go) %in% ge_observed_lowly_names))
names(ge_observed_lowly_genes) = names(all_go)


# ------------------------------------------------------------------------------------------------
#                         Molecular function GO Enrichment
# ------------------------------------------------------------------------------------------------

# nam physiological -----------------------------------------
np_observed_highly_go = new("topGOdata", allGenes = np_observed_highly_genes,
                            ontology="MF",
                            nodeSize = 10,
                            annot = annFUN.gene2GO,
                            gene2GO = all_go)

np_observed_lowly_go = new("topGOdata", allGenes = np_observed_lowly_genes,
                           ontology="MF",
                           nodeSize = 10,
                           annot = annFUN.gene2GO,
                           gene2GO = all_go)

np_observed_highly_enrich = runTest(np_observed_highly_go, algorithm = "classic", statistic="fisher")
np_observed_lowly_enrich = runTest(np_observed_lowly_go, algorithm = "classic", statistic="fisher")

# goodman expression ---------------------------------------
gp_observed_highly_go = new("topGOdata", allGenes = gp_observed_highly_genes,
                            ontology="MF",
                            nodeSize = 10,
                            annot = annFUN.gene2GO,
                            gene2GO = all_go)

gp_observed_lowly_go = new("topGOdata", allGenes = gp_observed_lowly_genes,
                           ontology="MF",
                           nodeSize = 10,
                           annot = annFUN.gene2GO,
                           gene2GO = all_go)

gp_observed_highly_enrich = runTest(gp_observed_highly_go, algorithm = "classic", statistic="fisher")
gp_observed_lowly_enrich = runTest(gp_observed_lowly_go, algorithm = "classic", statistic="fisher")

# Metabolite observed -------------------------------------
gm_observed_highly_go = new("topGOdata", allGenes = gm_observed_highly_genes,
                            ontology="MF",
                            nodeSize = 10,
                            annot = annFUN.gene2GO,
                            gene2GO = all_go)

gm_observed_lowly_go = new("topGOdata", allGenes = gm_observed_lowly_genes,
                           ontology="MF",
                           nodeSize = 10,
                           annot = annFUN.gene2GO,
                           gene2GO = all_go)

gm_observed_highly_enrich = runTest(gm_observed_highly_go, algorithm = "classic", statistic="fisher")
gm_observed_lowly_enrich = runTest(gm_observed_lowly_go, algorithm = "classic", statistic="fisher")

# expression observed ---------------------------------------
ge_observed_highly_go = new("topGOdata", allGenes = ge_observed_highly_genes,
                            ontology="MF",
                            nodeSize = 10,
                            annot = annFUN.gene2GO,
                            gene2GO = all_go)

ge_observed_lowly_go = new("topGOdata", allGenes = ge_observed_lowly_genes,
                           ontology="MF",
                           nodeSize = 10,
                           annot = annFUN.gene2GO,
                           gene2GO = all_go)

ge_observed_highly_enrich = runTest(ge_observed_highly_go, algorithm = "classic", statistic="fisher")
ge_observed_lowly_enrich = runTest(ge_observed_lowly_go, algorithm = "classic", statistic="fisher")

# aggregate results ------------------------------------------
mole_go_out = dplyr::bind_rows(
  GenTable(np_observed_highly_go, np_observed_highly_enrich, topNodes = length(score(np_observed_highly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "High") %>% 
    dplyr::mutate(Data = "NAM Physiological"),
  GenTable(np_observed_lowly_go, np_observed_lowly_enrich, topNodes = length(score(np_observed_lowly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "Low") %>% 
    dplyr::mutate(Data = "NAM Physiological"),
  GenTable(gp_observed_highly_go, gp_observed_highly_enrich, topNodes = length(score(gp_observed_highly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "High") %>% 
    dplyr::mutate(Data = "Goodman Physiological"),
  GenTable(gp_observed_lowly_go, gp_observed_lowly_enrich, topNodes = length(score(gp_observed_lowly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "Low") %>% 
    dplyr::mutate(Data = "Goodman Physiological"),
  GenTable(gm_observed_highly_go, gm_observed_highly_enrich, topNodes = length(score(gm_observed_highly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "High") %>% 
    dplyr::mutate(Data = "Goodman Metabolite"),
  GenTable(gm_observed_lowly_go, gm_observed_lowly_enrich, topNodes = length(score(gm_observed_lowly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "Low") %>% 
    dplyr::mutate(Data = "Goodman Metabolite"),
  GenTable(ge_observed_highly_go, ge_observed_highly_enrich, topNodes = length(score(ge_observed_highly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "High") %>% 
    dplyr::mutate(Data = "Goodman Expression"),
  GenTable(ge_observed_lowly_go, ge_observed_lowly_enrich, topNodes = length(score(ge_observed_lowly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "Low") %>% 
    dplyr::mutate(Data = "Goodman Expression")
)

# Format p-values
mole_go_out$result1 <- as.numeric(as.character(mole_go_out$result1))

# Add percentage of genes metric
mole_go_out$PercentageGenes <- (mole_go_out$Significant/mole_go_out$Annotated)*100
mole_go_out$PercentageGenesExpected <- (mole_go_out$Expected/mole_go_out$Annotated)*100

# Do FDR correction
mole_go_out$FDR <- p.adjust(mole_go_out$result1, method = "fdr")

# Add identifier on the type of GO enrichment this is
mole_go_out$go_type <- rep("MF", nrow(mole_go_out))


# --------------------------------------------------------------------------------------------------------------
#                       Biological GO Ontology
# --------------------------------------------------------------------------------------------------------------

# nam physiological ---------------------------------------------------------
np_observed_highly_go = new("topGOdata", allGenes = np_observed_highly_genes,
                            ontology="BP",
                            nodeSize = 10,
                            annot = annFUN.gene2GO,
                            gene2GO = all_go)

np_observed_lowly_go = new("topGOdata", allGenes = np_observed_lowly_genes,
                           ontology="BP",
                           nodeSize = 10,
                           annot = annFUN.gene2GO,
                           gene2GO = all_go)

np_observed_highly_enrich = runTest(np_observed_highly_go, algorithm = "classic", statistic="fisher")
np_observed_lowly_enrich = runTest(np_observed_lowly_go, algorithm = "classic", statistic="fisher")

# goodman expression --------------------------------------------------------
gp_observed_highly_go = new("topGOdata", allGenes = gp_observed_highly_genes,
                            ontology="BP",
                            nodeSize = 10,
                            annot = annFUN.gene2GO,
                            gene2GO = all_go)

gp_observed_lowly_go = new("topGOdata", allGenes = gp_observed_lowly_genes,
                           ontology="BP",
                           nodeSize = 10,
                           annot = annFUN.gene2GO,
                           gene2GO = all_go)

gp_observed_highly_enrich = runTest(gp_observed_highly_go, algorithm = "classic", statistic="fisher")
gp_observed_lowly_enrich = runTest(gp_observed_lowly_go, algorithm = "classic", statistic="fisher")

# Metabolite observed -------------------------------------------------------
gm_observed_highly_go = new("topGOdata", allGenes = gm_observed_highly_genes,
                            ontology="BP",
                            nodeSize = 10,
                            annot = annFUN.gene2GO,
                            gene2GO = all_go)

gm_observed_lowly_go = new("topGOdata", allGenes = gm_observed_lowly_genes,
                           ontology="BP",
                           nodeSize = 10,
                           annot = annFUN.gene2GO,
                           gene2GO = all_go)

gm_observed_highly_enrich = runTest(gm_observed_highly_go, algorithm = "classic", statistic="fisher")
gm_observed_lowly_enrich = runTest(gm_observed_lowly_go, algorithm = "classic", statistic="fisher")


# expression observed ------------------------------------------------------
ge_observed_highly_go = new("topGOdata", allGenes = ge_observed_highly_genes,
                            ontology="BP",
                            nodeSize = 10,
                            annot = annFUN.gene2GO,
                            gene2GO = all_go)

ge_observed_lowly_go = new("topGOdata", allGenes = ge_observed_lowly_genes,
                           ontology="BP",
                           nodeSize = 10,
                           annot = annFUN.gene2GO,
                           gene2GO = all_go)

ge_observed_highly_enrich = runTest(ge_observed_highly_go, algorithm = "classic", statistic="fisher")
ge_observed_lowly_enrich = runTest(ge_observed_lowly_go, algorithm = "classic", statistic="fisher")


# aggregate results ------------------------------------------------------
bio_go_out = dplyr::bind_rows(
  GenTable(np_observed_highly_go, np_observed_highly_enrich, topNodes = length(score(np_observed_highly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "High") %>% 
    dplyr::mutate(Data = "NAM Physiological"),
  GenTable(np_observed_lowly_go, np_observed_lowly_enrich, topNodes = length(score(np_observed_lowly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "Low") %>% 
    dplyr::mutate(Data = "NAM Physiological"),
  GenTable(gp_observed_highly_go, gp_observed_highly_enrich, topNodes = length(score(gp_observed_highly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "High") %>% 
    dplyr::mutate(Data = "Goodman Physiological"),
  GenTable(gp_observed_lowly_go, gp_observed_lowly_enrich, topNodes = length(score(gp_observed_lowly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "Low") %>% 
    dplyr::mutate(Data = "Goodman Physiological"),
  GenTable(gm_observed_highly_go, gm_observed_highly_enrich, topNodes = length(score(gm_observed_highly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "High") %>% 
    dplyr::mutate(Data = "Goodman Metabolite"),
  GenTable(gm_observed_lowly_go, gm_observed_lowly_enrich, topNodes = length(score(gm_observed_lowly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "Low") %>% 
    dplyr::mutate(Data = "Goodman Metabolite"),
  GenTable(ge_observed_highly_go, ge_observed_highly_enrich, topNodes = length(score(ge_observed_highly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "High") %>% 
    dplyr::mutate(Data = "Goodman Expression"),
  GenTable(ge_observed_lowly_go, ge_observed_lowly_enrich, topNodes = length(score(ge_observed_lowly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "Low") %>% 
    dplyr::mutate(Data = "Goodman Expression")
)

# Format p-values
bio_go_out$result1 <- as.numeric(as.character(bio_go_out$result1))

# Add percentage of genes metric
bio_go_out$PercentageGenes <- (bio_go_out$Significant/bio_go_out$Annotated)*100
bio_go_out$PercentageGenesExpected <- (bio_go_out$Expected/bio_go_out$Annotated)*100

# Do FDR correction
bio_go_out$FDR <- p.adjust(bio_go_out$result1, method = "fdr")

# Add identifier on the type of GO enrichment this is
bio_go_out$go_type <- rep("BP", nrow(bio_go_out))

# -----------------------------------------------------------------
# Cellular GO Ontology
# -----------------------------------------------------------------


# nam physiological ------------------------
np_observed_highly_go = new("topGOdata", allGenes = np_observed_highly_genes,
                            ontology="CC",
                            nodeSize = 10,
                            annot = annFUN.gene2GO,
                            gene2GO = all_go)

np_observed_lowly_go = new("topGOdata", allGenes = np_observed_lowly_genes,
                           ontology="CC",
                           nodeSize = 10,
                           annot = annFUN.gene2GO,
                           gene2GO = all_go)

np_observed_highly_enrich = runTest(np_observed_highly_go, algorithm = "classic", statistic="fisher")
np_observed_lowly_enrich = runTest(np_observed_lowly_go, algorithm = "classic", statistic="fisher")

# goodman expression -----------------------
gp_observed_highly_go = new("topGOdata", allGenes = gp_observed_highly_genes,
                            ontology="CC",
                            nodeSize = 10,
                            annot = annFUN.gene2GO,
                            gene2GO = all_go)

gp_observed_lowly_go = new("topGOdata", allGenes = gp_observed_lowly_genes,
                           ontology="CC",
                           nodeSize = 10,
                           annot = annFUN.gene2GO,
                           gene2GO = all_go)

gp_observed_highly_enrich = runTest(gp_observed_highly_go, algorithm = "classic", statistic="fisher")
gp_observed_lowly_enrich = runTest(gp_observed_lowly_go, algorithm = "classic", statistic="fisher")

# Metabolite observed ----------------------
gm_observed_highly_go = new("topGOdata", allGenes = gm_observed_highly_genes,
                            ontology="CC",
                            nodeSize = 10,
                            annot = annFUN.gene2GO,
                            gene2GO = all_go)

gm_observed_lowly_go = new("topGOdata", allGenes = gm_observed_lowly_genes,
                           ontology="CC",
                           nodeSize = 10,
                           annot = annFUN.gene2GO,
                           gene2GO = all_go)


gm_observed_highly_enrich = runTest(gm_observed_highly_go, algorithm = "classic", statistic="fisher")
gm_observed_lowly_enrich = runTest(gm_observed_lowly_go, algorithm = "classic", statistic="fisher")


# expression observed ----------------------
ge_observed_highly_go = new("topGOdata", allGenes = ge_observed_highly_genes,
                            ontology="CC",
                            nodeSize = 10,
                            annot = annFUN.gene2GO,
                            gene2GO = all_go)

ge_observed_lowly_go = new("topGOdata", allGenes = ge_observed_lowly_genes,
                           ontology="CC",
                           nodeSize = 10,
                           annot = annFUN.gene2GO,
                           gene2GO = all_go)

ge_observed_highly_enrich = runTest(ge_observed_highly_go, algorithm = "classic", statistic="fisher")
ge_observed_lowly_enrich = runTest(ge_observed_lowly_go, algorithm = "classic", statistic="fisher")


# aggregate results --------------
cell_go_out = dplyr::bind_rows(
  GenTable(np_observed_highly_go, np_observed_highly_enrich, topNodes = length(score(np_observed_highly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "High") %>% 
    dplyr::mutate(Data = "NAM Physiological"),
  GenTable(np_observed_lowly_go, np_observed_lowly_enrich, topNodes = length(score(np_observed_lowly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "Low") %>% 
    dplyr::mutate(Data = "NAM Physiological"),
  GenTable(gp_observed_highly_go, gp_observed_highly_enrich, topNodes = length(score(gp_observed_highly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "High") %>% 
    dplyr::mutate(Data = "Goodman Physiological"),
  GenTable(gp_observed_lowly_go, gp_observed_lowly_enrich, topNodes = length(score(gp_observed_lowly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "Low") %>% 
    dplyr::mutate(Data = "Goodman Physiological"),
  GenTable(gm_observed_highly_go, gm_observed_highly_enrich, topNodes = length(score(gm_observed_highly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "High") %>% 
    dplyr::mutate(Data = "Goodman Metabolite"),
  GenTable(gm_observed_lowly_go, gm_observed_lowly_enrich, topNodes = length(score(gm_observed_lowly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "Low") %>% 
    dplyr::mutate(Data = "Goodman Metabolite"),
  GenTable(ge_observed_highly_go, ge_observed_highly_enrich, topNodes = length(score(ge_observed_highly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "High") %>% 
    dplyr::mutate(Data = "Goodman Expression"),
  GenTable(ge_observed_lowly_go, ge_observed_lowly_enrich, topNodes = length(score(ge_observed_lowly_enrich))) %>%
    as_tibble() %>%
    dplyr::mutate(Effect = "Low") %>% 
    dplyr::mutate(Data = "Goodman Expression")
)

# Format p-values
cell_go_out$result1 <- as.numeric(as.character(cell_go_out$result1))

# Add percentage of genes metric
cell_go_out$PercentageGenes <- (cell_go_out$Significant/cell_go_out$Annotated)*100
cell_go_out$PercentageGenesExpected <- (cell_go_out$Expected/cell_go_out$Annotated)*100

# Do FDR correction
cell_go_out$FDR <- p.adjust(cell_go_out$result1, method = "fdr")

# Add identifier on the type of GO enrichment this is
cell_go_out$go_type <- rep("CC", nrow(cell_go_out))


# --------------------------------------------------------------------
#                  Set global plotting variables
# --------------------------------------------------------------------

axis_text_size <- 9.5
title_text_size <- 9.5
tag_size <- 8
annotate_size <- 1.5
pointSize <- 0.5
legend_key_size <- 0.8
legend_text <- 10
legend_text_size <- 10


# --------------------------------------------------------------------
#               Combine and Plot results
# --------------------------------------------------------------------

# check p-value cutoff
p_val_cutoff

# combine results from molecular function, biological process, and cellular components
go_results <- rbind(mole_go_out, bio_go_out, cell_go_out)

# Only gather significant go terms after FDR correction
go_results <- go_results[which(go_results$FDR <= p_val_cutoff), ]
dim(go_results)

# Combine names 
go_results$Term_Type <- paste0(go_results$Term, " (", go_results$go_type, ")")

# Correct names of populations
go_results$Data <- gsub("Physiological", "Field", go_results$Data)
go_results$Data <- gsub("Metabolite", "Mass Features", go_results$Data)

# Just gather terms in MF, BP, and CC separately and take the top 25 with the smallest corrected p-values
data <- tibble::as_tibble(go_results) %>%
  dplyr::group_by(go_type, Data) %>%
  dplyr::arrange(FDR, .by_group = TRUE) %>%
  dplyr::top_n(7)

# Plot results using ggplot2
# go_plot <- ggplot(data, aes(x = Effect, y = Term_Type)) + 
#   geom_point(size = 6, colour = "#999999") +
#   scale_size_area() +
#   facet_grid(.~factor(Data, levels=c("NAM Field","Goodman Field","Goodman Mass Features","Goodman Expression"))) +
#   labs(x="Pleiotropy Level", y="Gene Ontology Term") +
#   theme_bw() +
#   theme(axis.text=element_text(size=axis_text_size), 
#         axis.title = element_text(size=axis_text_size),
#         plot.title = element_text(size=title_text_size), 
#         legend.position="bottom", 
#         legend.text = element_text(size = legend_text_size),
#         legend.key.size = unit(legend_key_size, 'cm'),
#         strip.text = element_text(size = axis_text_size)) 

# Color points and change sizes
go_plot <- ggplot(data, aes(x = Effect, y = Term_Type)) + 
  geom_point(aes(size = PercentageGenes, colour = FDR)) +
  scale_size_area() +
  facet_grid(.~factor(Data, levels=c("NAM Field","Goodman Field","Goodman Mass Features","Goodman Expression"))) +
  labs(x="Pleiotropy Level", y="Gene Ontology Term") +
  theme_bw() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.key.size = unit(legend_key_size, 'cm'),
        strip.text = element_text(size = axis_text_size)) 

go_plot  

# Save plot
ggsave("~/git_projects/pleiotropy/images/go_highly_lowly_top_bottom_25_percentile_adjusted_pleiotropy_above0.png",
       plot = go_plot,
       width = 10,
       height = 6,
       units = "in")
