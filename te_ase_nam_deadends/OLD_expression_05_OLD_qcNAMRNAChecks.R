# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-12-21
#
# Description 
#   - 07_nam_rna_checks.R
#   - Normalize feature count data from HTSeq
# Detailed Purpose:
#    The main purpose of this Rscript is to generate QC figures for
#    RNA-seq alignment quality of forward reads from the NAM
#    h RNybridA-seq experiment. Count data is read from the result
#    of the prior script (`06_count_normalization.R`).
#--------------------------------------------------------------------

# === Preamble ======================================================

## Load packages ----
library(data.table)
library(dplyr)
library(ggplot2)
library(magrittr)
library(readr)
library(tidyr)


# === QC - Estimate saturation ======================================

### Count data
count_raw <- data.table::fread("/workdir/mbb262/nam_hybrid_rnaseq/output/counts/nam_counts_raw.csv")

## Convert count data to matrix ----
count_mat <- count_raw[, -1] %>% as.matrix()
rownames(count_mat) <- count_raw[,1] %>% as.matrix()

## Use RNAseQC package for saturation curves ----
# Note: this package does not compile so I took the functions I needed directly
#       from github and source them here
source("~/git_projects/te_ase_nam/src/07_estimate_rna_saturation_plot.R")

### NOTE: this can be time-consuming...
sat_path <- "/workdir/mbb262/nam_hybrid_rnaseq/output/counts/nam_saturation_values.csv"
if (!file.exists(sat_path)) {
    ndepth <- 15
    sat <- estimate_saturation(
        counts = count_mat,
        min_counts = 2,
        verbose = TRUE,
        ndepths = ndepth
    )
    write.csv(
        x = sat,
        file = sat_path,
        row.names = FALSE
    )
} else {
    sat <- readr::read_csv(file = sat_path)
}


# === QC - Visualize data (saturation) ==============================

## Visualize saturation ----
satplot <- sat %>%
    ggplot() +
    aes(x = depth, y = sat, color = sample) +
    labs(x = "Reads", y = "Genes above threshold") +
    geom_line(size = 1.2, alpha = 0.9) +
    xlab("Number of reads") +
    ylab("# of genes with at least 2 reads") +
    scale_x_continuous(labels = function(l) {
        paste0(round(l/1e6,1),"M")
    }) +
    theme(legend.position = "none", legend.title = element_blank())
satplot


# === QC - Visualize data (hisat2 alignments) =========================

# Load in data
align <- read.delim("~/Box Sync/Cornell_PhD/labProjects/nam_hybrid_te/data/multiqc/multiqc_hisat2.txt", sep = "\t")
align <- read.delim("/workdir/mbb262/multiqc_data/multiqc_hisat2.txt", sep = "\t")

# Pick samples that map less than 50% overall
align %>% 
    select("Sample", "overall_alignment_rate") %>% 
    filter(overall_alignment_rate <= 50)

## Convert data to long format ----
long_hisat_data <- align %>%
    dplyr::arrange(dplyr::desc(paired_total)) %>%
    dplyr::mutate(sample_id = seq_len(nrow(align) %>% as.numeric())) %>%
    tidyr::gather(key = "map_type", value = "value", c(3:6, 8:10)) %>%
    dplyr::mutate(
        map_type = factor(
            x = map_type,
            levels = c(
                "paired_aligned_none",
                "paired_aligned_one", 
                "paired_aligned_multi",
                "paired_aligned_discord_one",
                "unpaired_aligned_none",
                "unpaired_aligned_one",
                "unpaired_aligned_multi"
            )
        )
    )


## Plot hisat alignment scores ----
hisat_scores <- long_hisat_data %>%
    ggplot() +
    aes(x = sample_id, y = value, fill = map_type) +
    geom_area(position = position_stack(reverse = TRUE)) +
    scale_fill_manual(
        values = c("#002240", "#4579A6", "#77B4EE", "#BA0948", "#E1AB76", "#E3D378", "#780505"),
        name = "",
        breaks = c(
            "paired_aligned_none",
            "paired_aligned_one",
            "paired_aligned_multi",
            "paired_aligned_discord_one",
            "unpaired_aligned_none",
            "unpaired_aligned_one",
            "unpaired_aligned_multi"
        ),
        labels = c(
            "PE, aligned 0 times",
            "PE, aligned exactly 1 time",
            "PE, aligned >1 times",
            "PE, aligned discordantly 1 time",
            "Unpaired, aligned exactly 1 time",
            "Unpaired, aligned >1 times",
            "Unpaired, aligned 0 times"
        )
    ) +
    xlab("NAM inbreds and hybrids arranged by total input reads") +
    ylab("Number of reads") +
    ggtitle("HISAT2 alignment scores") +
    expand_limits( x = c(0,NA), y = c(0,NA)) +
    scale_y_continuous(labels = function(l) {
        paste0(round(l / 1e6, 1), "M")
    }) +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    )


# === Write plots to disk ============================================


## Saturation curve ----
ggsave(
    filename = "~/git_projects/te_ase_nam/figs/saturation.png",
    plot = satplot,device = "png",width = 8, height = 5, dpi = "retina"
)


## hisat alignments ----
ggsave(
    filename = "~/git_projects/te_ase_nam/figs/hisat_scores.png",
    plot = hisat_scores, device = "png",width = 10, height = 5, dpi = "retina"
)


# === Run PCA to check sample placement =============================

## DEseq object is from script 06
pcaData <- DESeq2::plotPCA(ddsVST, intgroup=c("cultivar", "tissue", "age"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Make plot, age by tissue
age_tissue <- ggplot(pcaData, aes(PC1, PC2, color=age, shape=tissue)) +
    geom_point(size=2) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()

# Make plot, age by tissue
ggplot(pcaData, aes(PC1, PC2, color=age, shape=tissue, label = cultivar)) +
    geom_point(size=2) +
    geom_text(hjust=0, vjust=0, size = 3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()


# === Write plots to disk ============================================


## PCA: age by tissue ----
ggsave(
    filename = "~/git_projects/te_ase_nam/figs/pca_age_tissue.png",
    plot = age_tissue,device = "png",width = 8, height = 5, dpi = "retina"
)


# # === Filter bad stuff and visualize again (hisat stats) =============================
# 
# ## Filter alignment data ----
# hisat_data_filt <- align %>%
#     dplyr::filter(paired_total >= 1e6) %>%
#     dplyr::filter(unmapped_short / paired_total <= 0.5)
# 
# 
# ## Convert data to long format ----
# long_hisat_data_filt <- hisat_data_filt %>%
#     dplyr::arrange(dplyr::desc(paired_total)) %>%
#     dplyr::mutate(sample_id = seq_len(nrow(hisat_data_filt) %>% as.numeric())) %>%
#     tidyr::gather(key = "map_type", value = "value", 4:8) %>%
#     dplyr::mutate(
#         map_type = factor(
#             x = map_type,
#             levels = c(
#                 "paired_aligned_none",
#                 "paired_aligned_one", 
#                 "paired_aligned_multi",
#                 "paired_aligned_discord_one",
#                 "unpaired_aligned_none",
#                 "unpaired_aligned_one",
#                 "unpaired_aligned_multi"
#             )
#         )
#     )
# 
# ## Plot hisat alignment scores ----
# hisat_scores_filt <- long_hisat_data_filt %>%
#     ggplot() +
#     aes(x = sample_id, y = value, fill = map_type) +
#     geom_area(position = position_stack(reverse = TRUE)) +
#     scale_fill_manual(
#         values = c("#002240", "#4579A6", "#77B4EE", "#BA0948", "#E1AB76", "#E3D378", "#780505"),
#         name = "",
#         breaks = c(
#             "paired_aligned_none",
#             "paired_aligned_one",
#             "paired_aligned_multi",
#             "paired_aligned_discord_one",
#             "unpaired_aligned_none",
#             "unpaired_aligned_one",
#             "unpaired_aligned_multi"
#         ),
#         labels = c(
#             "PE, aligned 0 times",
#             "PE, aligned exactly 1 time",
#             "PE, aligned >1 times",
#             "PE, aligned discordantly 1 time",
#             "Unpaired, aligned exactly 1 time",
#             "Unpaired, aligned >1 times",
#             "Unpaired, aligned 0 times"
#         )
#     ) +
#     xlab("NAM inbreds and hybrids arranged by total input reads") +
#     ylab("Number of reads") +
#     ggtitle("HISAT2 alignment scores (filtered)") +
#     expand_limits( x = c(0,NA), y = c(0,NA)) +
#     scale_y_continuous(labels = function(l) {
#         paste0(round(l / 1e6, 1), "M")
#     }) +
#     theme_minimal() +
#     theme(
#         legend.position = "bottom",
#         axis.text.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank()
#     )
# 
# 
# ## Which taxa IDs have poor quality? ----
# 
# ### NOTE: This data will be used in script `06`
# bad_taxa <- dplyr::setdiff(
#     x = hisat_data$sample,
#     y = hisat_data_filt$sample
# )
# 
# 
# ## Save ----
# 
# ### Plot
# ggsave(
#     filename = "figs/hisat_scores_filtered.png",
#     plot = hisat_scores_filt,
#     device = "png",
#     width = 10,
#     height = 5,
#     dpi = "retina"
# )
# 
# ### Bad taxa IDs
# readr::write_lines(x = bad_taxa, path = "local_data/bad_taxa.txt")
# 

