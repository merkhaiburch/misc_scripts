# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-08-29
#
# Description 
#   - Exploring maize TE annotations from Stitzer 2019
# ---------------------------------------------------------------

# packages
library(dplyr)
library(ggplot2)
library(patchwork)
library(data.table)

# Load data
te1 <- data.table::fread("git_projects/te_regulation/data/stitzer_figshare/B73.LTRAGE.allDescriptors.2019-01-31.txt")
closest_gene <- data.table::fread("git_projects/te_regulation/data/stitzer_figshare/B73_closest_gene.2019-10-21.txt")

# Intersect-Join te info with gene information
data.table::setkey(te1, TEID)
data.table::setkey(closest_gene, TEID)
closest_te1 <- te1[closest_gene, nomatch = 0]

# How many familes are there?
length(unique(te1$fam))

# ---------------------------------
# What TEs are within 3kb of a TSS?
# ---------------------------------

# Subset TEs to those within 3kb and not within a gene and 
# A bunch of familes are only seen once, try removing those, counting how many familes are left
closest_te1_3kb <- closest_te1 %>% 
  filter(closest <= 3000, closest >= 1) %>% 
  group_by(fam) %>% 
  filter(n() > 1)  # require each level with more than 4 obs

# count number of unique families
length(unique(closest_te1_3kb$fam)) # 846
length(unique(te1$fam)) # Total = 27k

# Subset file to only gene expression columns
closest_te1_3kb <- closest_te1_3kb %>% 
  select_if(!grepl("TEfam_", names(.)) & !grepl("_flank_", names(.)))

# Look at gene_Mature_Leaf_8
plot(log10(closest_te1_3kb$gene_Mature_Leaf_8+1), closest_te1_3kb$age)

# Age, model is sig
model_age <- lm(log10(gene_Mature_Leaf_8+1) ~ age, data = closest_te1_3kb)
summary(model_age)

# age + superfamily, model sig
model_age_sup <- lm(log10(gene_Mature_Leaf_8+1) ~ age + sup, data = closest_te1_3kb)
summary(model_age_sup)

# age + family, model ns
model_age_fam <- lm(log10(gene_Mature_Leaf_8+1) ~ age + fam, data = closest_te1_3kb)
summary(model_age_fam)

# Count TE families in this area, is there an enrichment? Yes, most only represented once
count_families <- data.frame(table(closest_te1_3kb$fam))
hist(count_families$Freq)


# ------------------------------
# Plotting afe+fam+sup+distance
# ------------------------------

# All TEs genome wide
# closest gene + superfamily, model sig
temp1 <- closest_te1 %>% filter(closest >= 1)
model_dist_sup <- lm(log10(gene_Mature_Leaf_8+1) ~ closest+ sup, data = temp1)
summary(model_dist_sup)
a <- ggplot(temp1, aes(x = closest, y = log10(gene_Mature_Leaf_8+1), color = sup))+
  geom_point()+
  geom_smooth(method = "lm")  +
  ylab("log10(B73 Mature Leaf Expression)") +
  xlab("Closest Gene") +
  ggtitle("B73 Leaf Expression ~ Distance from gene + \n all TE superfamilies")


# TEs within 3kb of a gene
temp <- closest_te1 %>%
  filter(closest <= 3000, closest >= 1)
length(unique(temp$fam))

# distance + superfamily, model sig
model_dist_sup <- lm(log10(gene_Mature_Leaf_8+1) ~ closest+  sup, data = temp)
summary(model_dist_sup)

b <- ggplot(temp, aes(x = closest, y = log10(gene_Mature_Leaf_8+1), color = sup))+
  geom_point()+
  geom_smooth(method = "lm") +
  ylab("log10(B73 Mature Leaf Expression)") +
  xlab("Closest Gene") +
  ggtitle("B73 Leaf Expression ~ Distance from gene + \n TE superfamilies <3Kb to a gene")


# Plot together
a+b

# # Facet out by superfamily
# c <- ggplot(te1, aes(x = closest, y = log10(gene_Mature_Leaf_8+1), color = sup))+
#   geom_point()+
#   geom_smooth(method = "lm") +
#   facet_wrap(~sup, ncol = 1) +
#   d <- ggplot(temp, aes(x = closest, y = log10(gene_Mature_Leaf_8+1), color = sup))+
#   geom_point()+
#   geom_smooth(method = "lm") +
#   facet_wrap(~sup, ncol = 1)
# c+d


# Density plots of number of TE familes away from genes
temp0 <- te1 %>% select(TEID, fam, sup, closest) %>% filter(closest <= 500000, closest >= 1) %>% group_by(fam) %>% mutate(count = n()) %>% distinct(fam, .keep_all = TRUE)
temp1 <- te1 %>% select(TEID, fam, sup, closest) %>% filter(closest <= 100000, closest >= 1) %>% group_by(fam) %>% mutate(count = n()) %>% distinct(fam, .keep_all = TRUE)
temp2 <- te1 %>% select(TEID, fam, sup, closest) %>% filter(closest <= 50000, closest >= 1) %>% group_by(fam) %>% mutate(count = n()) %>% distinct(fam, .keep_all = TRUE)
temp3 <- te1 %>% select(TEID, fam, sup, closest) %>% filter(closest <= 30000, closest >= 1) %>% group_by(fam) %>% mutate(count = n()) %>% distinct(fam, .keep_all = TRUE)
temp4 <- te1 %>% select(TEID, fam, sup, closest) %>% filter(closest <= 20000, closest >= 1) %>% group_by(fam) %>% mutate(count = n()) %>% distinct(fam, .keep_all = TRUE)
temp5 <- te1 %>% select(TEID, fam, sup, closest) %>% filter(closest <= 10000, closest >= 1) %>% group_by(fam) %>% mutate(count = n()) %>% distinct(fam, .keep_all = TRUE)
temp6 <- te1 %>% select(TEID, fam, sup, closest) %>% filter(closest <= 3000, closest >= 1) %>% group_by(fam) %>% mutate(count = n()) %>% distinct(fam, .keep_all = TRUE)
length(unique(temp6$fam))
# add names
temp0$distance <- rep(500000, nrow(temp0))
temp1$distance <- rep(100000, nrow(temp1))
temp2$distance <- rep(50000, nrow(temp2))
temp3$distance <- rep(30000, nrow(temp3))
temp4$distance <- rep(20000, nrow(temp4))
temp5$distance <- rep(10000, nrow(temp5))
temp6$distance <- rep(3000, nrow(temp6))

# Rbind all together
together <- rbind(temp0, temp1, temp2, temp3, temp4, temp5, temp6) %>% filter(count >20)


# -------------
# Analysis
#--------------

# Look at relationship between B73 mature leaf expression and 
#     all TEs genome wide + their distance away from genes
model_dist_sup <- lm(log10(gene_Mature_Leaf_8+1) ~ closest+ sup, data = te1)
summary(model_dist_sup)

# Plot this model
ggplot(te1, aes(x = closest, y = log10(gene_Mature_Leaf_8+1), color = sup))+
  geom_point()+
  geom_smooth(method = "lm")  +
  ylab("log10(B73 Mature Leaf Expression)") +
  ggtitle("B73 Leaf Expression ~ Distance from gene + all TE superfamilies")

