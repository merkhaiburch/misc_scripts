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
sweeps <- data.table::fread("../Downloads/sweeps.csv", header = TRUE)
colnames(sweeps) <- c("seqnames", "start", "end", "xpclr", "comparison")

# Load in pleiotropic entropy data
entropy <- data.table::fread("../Downloads/goodman_nam_allTraitTypes_entropy_gwaCounts.txt")

# interval data
ranges <- read.csv("../git_projects/pleiotropy/data/genic_intergenic_intervals_b73_v4.49.csv")

# Combine entropy intervals with coordinates
entropy <- merge(x = entropy, y = ranges, by.x= "Interval", by.y = "rr_id") |> data.table()

# find overlaps between sweeps and my intervals
setkey(sweeps, seqnames, start, end)
entropy <- foverlaps( entropy, sweeps, by.x = c("seqnames", "start", "end"), mult = "first")
dim(entropy)
length(unique(entropy$Interval))

# Join back with ranges
# entropy$xpclr[is.na(entropy$xpclr)] <- 0

# rearrange
entropy <- entropy[,c(6:12, 1:5, 13:14)] |> data.frame()


# ------------------------
# Create bins for entrpoy
# ------------------------

# shannon_entropy_nam_filtered_all
q <- quantile(entropy[,2], seq(0,1, .333))
entropy$binned_nfa <- entropy[,2]
entropy$binned_nfa[entropy$binned_nfa >= q[1] & entropy$binned_nfa <= q[2]] <-1
entropy$binned_nfa[entropy$binned_nfa > q[2] & entropy$binned_nfa <= q[3]] <- 2
entropy$binned_nfa[entropy$binned_nfa > q[3] & entropy$binned_nfa <= max(entropy[,2])] <- 3

#shannon_entropy_nam_random_all
q <- quantile(entropy[,3], seq(0,1, .333))
entropy$binned_nra <- entropy[,3]
entropy$binned_nra[entropy$binned_nra >= q[1] & entropy$binned_nra <= q[2]] <- 1
entropy$binned_nra[entropy$binned_nra > q[2] & entropy$binned_nra <= q[3]] <- 2
entropy$binned_nra[entropy$binned_nra > q[3] & entropy$binned_nra <= max(entropy[,3])] <- 3

#shannon_entropy_goodman_filtered_metabolite
q <- quantile(entropy[,4], seq(0,1, .333))
entropy$binned_gfm <- entropy[,4]
entropy$binned_gfm[entropy$binned_gfm >= q[1] & entropy$binned_gfm <= q[2]] <- 1
entropy$binned_gfm[entropy$binned_gfm > q[2] & entropy$binned_gfm <= q[3]] <- 2
entropy$binned_gfm[entropy$binned_gfm > q[3] & entropy$binned_gfm <= max(entropy[,4])] <- 3

#shannon_entropy_goodman_random_metabolite
q <- quantile(entropy[,5], seq(0,1, .333))
entropy$binned_grm <- entropy[,5]
entropy$binned_grm[entropy$binned_grm >= q[1] & entropy$binned_grm <= q[2]] <- 1
entropy$binned_grm[entropy$binned_grm > q[2] & entropy$binned_grm <= q[3]] <- 2
entropy$binned_grm[entropy$binned_grm > q[3] & entropy$binned_grm <= max(entropy[,5])] <- 3

# shannon_entropy_goodman_filtered_physiological
q <- quantile(entropy[,6], seq(0,1, .333))
entropy$binned_gfp <- entropy[,6]
entropy$binned_gfp[entropy$binned_gfp >= q[1] & entropy$binned_gfp <= q[2]] <- 1
entropy$binned_gfp[entropy$binned_gfp > q[2] & entropy$binned_gfp <= q[3]] <- 2
entropy$binned_gfp[entropy$binned_gfp > q[3] & entropy$binned_gfp <= max(entropy[,6])] <- 3

# shannon_entropy_goodman_random_physiological
q <- quantile(entropy[,7], seq(0,1, .333))
entropy$binned_grp <- entropy[,7]
entropy$binned_grp[entropy$binned_grp >= q[1] & entropy$binned_grp <= q[2]] <- 1
entropy$binned_grp[entropy$binned_grp > q[2] & entropy$binned_grp <= q[3]] <- 2
entropy$binned_grp[entropy$binned_grp > q[3] & entropy$binned_grp <= max(entropy[,7])] <- 3


# melt dataframe
melted_entropy <- entropy |> tidyr::pivot_longer(
  cols = binned_nfa:binned_grp,
  names_prefix = "binned_",
  values_to = "binned_level")
# 
melted_entropy <- melted_entropy |> tidyr::pivot_longer(
  cols = shannon_entropy_nam_filtered_all:shannon_entropy_goodman_random_physiological,
  names_to = "pop_trait",
  names_prefix = "shannon_entropy_",
  values_to = "entropy")

# Plot
ggplot(melted_entropy, aes(x=as.factor(binned_level), y=xpclr, fill=comparison)) +
  geom_boxplot() +
  facet_grid(cols = vars(pop_trait))


# Subset out bu pop-trait
# nam all
temp_nam <-  tidyr::pivot_longer(melted_entropy[,c(1,2:3, 11:12)],
  cols = shannon_entropy_nam_filtered_all:shannon_entropy_nam_random_all,
  names_to = "pop_trait",
  names_prefix = "shannon_entropy_",
  values_to = "entropy")

ggplot(temp_nam, aes(x=xpclr, y=entropy, colour = comparison)) +
  geom_smooth(method = "lm", aes(linetype = pop_trait))


# metabolites
temp_meta <-  tidyr::pivot_longer(melted_entropy[,c(1,4:5, 11:12)],
                             cols = shannon_entropy_goodman_filtered_metabolite:shannon_entropy_goodman_random_metabolite,
                             names_to = "pop_trait",
                             names_prefix = "shannon_entropy_",
                             values_to = "entropy")

ggplot(temp_meta, aes(x=xpclr, y=entropy, colour = comparison)) +
  geom_smooth(method = "lm", aes(linetype = pop_trait))


# goodman phys
temp_gp <-  tidyr::pivot_longer(melted_entropy[,c(1,6:7, 11:12)],
                             cols = shannon_entropy_goodman_filtered_physiological:shannon_entropy_goodman_random_physiological,
                             names_to = "pop_trait",
                             names_prefix = "shannon_entropy_",
                             values_to = "entropy")

ggplot(temp_gp, aes(x=xpclr, y=entropy, colour = comparison)) +
  geom_smooth(method = "lm", aes(linetype = pop_trait))
