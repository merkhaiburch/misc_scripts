# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-10-06 
#
# Description 
#   - Querying if there is an enrichment of sweep sites in different
#   - levels/classifications of pleiotropic sites
# ---------------------------------------------------------------

# Load in packages
library(magrittr)
library(tidyr)
library(data.table)
library(ggplot2)


# -------------------------------
# Data gathering & formatting
# -------------------------------

# Load in sweep data from:
# Evolutionary and functional genomics of DNA methylation in maize domestication and improvement. 
sweeps <- data.table::fread("/workdir/mbb262/interval_data/sweeps.csv", header = TRUE)
colnames(sweeps) <- c("seqid", "start", "end", "xpclr", "comparison")

# Load in pleiotropic entropy data
all_interval_data <- data.table::fread("/workdir/mbb262/interval_data/all_interval_data.csv")
all_interval_data_melted <- data.table::fread("/workdir/mbb262/interval_data/all_interval_data_melted.csv")

# find overlaps between sweeps and my intervals
setkey(sweeps, seqid, start, end)
entropy <- foverlaps(all_interval_data_melted, sweeps, by.x = c("seqid", "start", "end"), mult = "first")
dim(entropy)
length(unique(entropy$rr_id))

# entropy$xpclr <- log(entropy$xpclr)

# ------------------------------------
#   Set global plotting variables
# ------------------------------------

axis_text_size <- 5
title_text_size <- 6.5
tag_size <- 9
filt_text <- 11
ran_text <- 9
r2_text_x <- 45
annotate_size <- 1.5
pointSize <- 0.5
legend_key_size <- 0.3
legend_text <- 2.5
lineSize <- 0.8


# nam physiological-----------------------------------------------------------------------------------

filtered_df <- entropy %>% filter(pop_trait == "NAM Physiological") %>% 
  select(rr_id, range_type, fitted_entropy, xpclr, filtering_level, pop_trait, comparison)

xpclr1 <- ggplot(filtered_df, aes(x = xpclr, y = fitted_entropy, color = comparison)) +
  geom_point(alpha = .5, size = pointSize) +
  geom_smooth(method = lm, size = pointSize, aes(linetype = filtering_level)) + 
  labs(x = "XPCLR", 
       y = "Adjusted Entropy", 
       title = "NAM Physiological") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text = element_text(size = legend_text),
        legend.key.size = unit(legend_key_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  labs(tag = "A") +
  ylim(-12,12)


# goodman physiological--------------------------------------------------------------------------------

filtered_df <- entropy %>% 
  filter(pop_trait == "Goodman Physiological")%>% 
  select(rr_id, range_type, fitted_entropy, xpclr, filtering_level, comparison)

xpclr2 <- ggplot(filtered_df, aes(x = xpclr, y = fitted_entropy, color = comparison)) +
  geom_point(alpha = .5, size = pointSize) +
  geom_smooth(method = lm, size = pointSize, aes(linetype = filtering_level)) + 
  labs(x = "XPCLR", 
       y = "Adjusted Entropy", 
       title = "Goodman Physiological") +
  theme_classic()  +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text = element_text(size = legend_text),
        legend.key.size = unit(legend_key_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  labs(tag = "B") +
  ylim(-12,12)


# goodman metabolite---------------------------------------------------------------------------------

filtered_df <- entropy %>% 
  filter(pop_trait == "Goodman Metabolite") %>% 
  select(rr_id, range_type, fitted_entropy, xpclr, filtering_level, comparison)

xpclr3 <- ggplot(filtered_df, aes(x = xpclr, y = fitted_entropy, color = comparison)) +
  geom_point(alpha = .5, size = pointSize) +
  geom_smooth(method = lm, size = pointSize, aes(linetype = filtering_level)) + 
  labs(x = "XPCLR", 
       y = "Adjusted Entropy", 
       title = "Goodman Metabolite") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text = element_text(size = legend_text),
        legend.key.size = unit(legend_key_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  labs(tag = "C") +
  ylim(-12,12)


# goodman expression---------------------------------------------------------------------------------

filtered_df <- entropy %>% 
  filter(pop_trait == "Goodman Expression") %>% 
  select(rr_id, range_type, fitted_entropy, xpclr, filtering_level, comparison)

xpclr4 <- ggplot(filtered_df, aes(x = xpclr, y = fitted_entropy, color = comparison)) +
  geom_point(alpha = .5, size = pointSize) +
  geom_smooth(method = lm, size = pointSize, aes(linetype = filtering_level)) + 
  labs(x = "XPCLR", 
       y = "Adjusted Entropy", 
       title = "Goodman Expression") +
  theme_classic()  +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text = element_text(size = legend_text),
        legend.key.size = unit(legend_key_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  labs(tag = "D") +
  ylim(-12,12)


# Test plot -------------------------------------------------------------------------------------------
# (xpclr1 + xpclr2)/(xpclr3 + xpclr4)


# export plot ----------------------------------------------------------------------------------------

# ggsave("/home/mbb262/git_projects/pleiotropy/images/adjusted_entropy_vs_xpclr.tiff",
#        plot = (xpclr1 + xpclr2)/(xpclr3 + xpclr4),
#        width = 6.5,
#        height = 4,
#        units = "in")

ggsave("/home/mbb262/git_projects/pleiotropy/images/adjusted_entropy_vs_xpclr.png",
       plot = (xpclr1 + xpclr2)/(xpclr3 + xpclr4),
       width = 7,
       height = 4, 
       units = "in")


# --------------------------------
# pullout expression data, bin it, see if it has an effect

temp <- entropy %>% 
  filter(pop_trait == "Goodman Expression") %>% 
  select(rr_id, range_type, fitted_entropy, xpclr, filtering_level, comparison)
temp1 <- temp %>% filter(fitted_entropy >= 2)
dim(temp1)

a<- ggplot(temp1, aes(x=xpclr, y = fitted_entropy, color=comparison)) +
  geom_point(alpha = .5) +
  geom_smooth(position="dodge", method = "lm", aes(linetype = filtering_level))
a

# 
# 
# 
# # ------------------------
# # Create bins for entrpoy
# # ------------------------
# 
# # shannon_entropy_nam_filtered_all
# q <- quantile(entropy[,2], seq(0,1, .333))
# entropy$binned_nfa <- entropy[,2]
# entropy$binned_nfa[entropy$binned_nfa >= q[1] & entropy$binned_nfa <= q[2]] <-1
# entropy$binned_nfa[entropy$binned_nfa > q[2] & entropy$binned_nfa <= q[3]] <- 2
# entropy$binned_nfa[entropy$binned_nfa > q[3] & entropy$binned_nfa <= max(entropy[,2])] <- 3
# 
# #shannon_entropy_nam_random_all
# q <- quantile(entropy[,3], seq(0,1, .333))
# entropy$binned_nra <- entropy[,3]
# entropy$binned_nra[entropy$binned_nra >= q[1] & entropy$binned_nra <= q[2]] <- 1
# entropy$binned_nra[entropy$binned_nra > q[2] & entropy$binned_nra <= q[3]] <- 2
# entropy$binned_nra[entropy$binned_nra > q[3] & entropy$binned_nra <= max(entropy[,3])] <- 3
# 
# #shannon_entropy_goodman_filtered_metabolite
# q <- quantile(entropy[,4], seq(0,1, .333))
# entropy$binned_gfm <- entropy[,4]
# entropy$binned_gfm[entropy$binned_gfm >= q[1] & entropy$binned_gfm <= q[2]] <- 1
# entropy$binned_gfm[entropy$binned_gfm > q[2] & entropy$binned_gfm <= q[3]] <- 2
# entropy$binned_gfm[entropy$binned_gfm > q[3] & entropy$binned_gfm <= max(entropy[,4])] <- 3
# 
# #shannon_entropy_goodman_random_metabolite
# q <- quantile(entropy[,5], seq(0,1, .333))
# entropy$binned_grm <- entropy[,5]
# entropy$binned_grm[entropy$binned_grm >= q[1] & entropy$binned_grm <= q[2]] <- 1
# entropy$binned_grm[entropy$binned_grm > q[2] & entropy$binned_grm <= q[3]] <- 2
# entropy$binned_grm[entropy$binned_grm > q[3] & entropy$binned_grm <= max(entropy[,5])] <- 3
# 
# # shannon_entropy_goodman_filtered_physiological
# q <- quantile(entropy[,6], seq(0,1, .333))
# entropy$binned_gfp <- entropy[,6]
# entropy$binned_gfp[entropy$binned_gfp >= q[1] & entropy$binned_gfp <= q[2]] <- 1
# entropy$binned_gfp[entropy$binned_gfp > q[2] & entropy$binned_gfp <= q[3]] <- 2
# entropy$binned_gfp[entropy$binned_gfp > q[3] & entropy$binned_gfp <= max(entropy[,6])] <- 3
# 
# # shannon_entropy_goodman_random_physiological
# q <- quantile(entropy[,7], seq(0,1, .333))
# entropy$binned_grp <- entropy[,7]
# entropy$binned_grp[entropy$binned_grp >= q[1] & entropy$binned_grp <= q[2]] <- 1
# entropy$binned_grp[entropy$binned_grp > q[2] & entropy$binned_grp <= q[3]] <- 2
# entropy$binned_grp[entropy$binned_grp > q[3] & entropy$binned_grp <= max(entropy[,7])] <- 3
# 
# 
# # melt dataframe
# melted_entropy <- entropy |> tidyr::pivot_longer(
#   cols = binned_nfa:binned_grp,
#   names_prefix = "binned_",
#   values_to = "binned_level")
# # 
# melted_entropy <- melted_entropy |> tidyr::pivot_longer(
#   cols = shannon_entropy_nam_filtered_all:shannon_entropy_goodman_random_physiological,
#   names_to = "pop_trait",
#   names_prefix = "shannon_entropy_",
#   values_to = "entropy")
# 
# # Plot
# ggplot(melted_entropy, aes(x=as.factor(binned_level), y=xpclr, fill=comparison)) +
#   geom_boxplot() +
#   facet_grid(cols = vars(pop_trait))
# 
# 
# # Subset out bu pop-trait
# # nam all
# temp_nam <-  tidyr::pivot_longer(melted_entropy[,c(1,2:3, 11:12)],
#   cols = shannon_entropy_nam_filtered_all:shannon_entropy_nam_random_all,
#   names_to = "pop_trait",
#   names_prefix = "shannon_entropy_",
#   values_to = "entropy")
# 
# ggplot(temp_nam, aes(x=xpclr, y=entropy, colour = comparison)) +
#   geom_smooth(method = "lm", aes(linetype = pop_trait))
# 
# 
# # metabolites
# temp_meta <-  tidyr::pivot_longer(melted_entropy[,c(1,4:5, 11:12)],
#                              cols = shannon_entropy_goodman_filtered_metabolite:shannon_entropy_goodman_random_metabolite,
#                              names_to = "pop_trait",
#                              names_prefix = "shannon_entropy_",
#                              values_to = "entropy")
# 
# ggplot(temp_meta, aes(x=xpclr, y=entropy, colour = comparison)) +
#   geom_smooth(method = "lm", aes(linetype = pop_trait))
# 
# 
# # goodman phys
# temp_gp <-  tidyr::pivot_longer(melted_entropy[,c(1,6:7, 11:12)],
#                              cols = shannon_entropy_goodman_filtered_physiological:shannon_entropy_goodman_random_physiological,
#                              names_to = "pop_trait",
#                              names_prefix = "shannon_entropy_",
#                              values_to = "entropy")
# 
# ggplot(temp_gp, aes(x=xpclr, y=entropy, colour = comparison)) +
#   geom_smooth(method = "lm", aes(linetype = pop_trait))
