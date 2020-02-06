# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-02-03 
#
# Description 
#   - Test effect of collapsing PHG haplotypes with 26 NAM founders
#   - maxDiv levels = 0.0001, 0.0005, 0.005, 0.001, 0.01
# ---------------------------------------------------------------

# Set workdir
setwd("~/Box Sync/Cornell_PhD/labProjects/hackathons/2020_02_03_hackathon/data")

# Libraries
library(ggplot2)
library(dplyr)


# ------------------
# Evan's dataframe
# ------------------

# Taxa ID by haplotype seqeunce length
ggplot(hap, aes(x = as.factor(gamete_grp_id), y = seq_len)) + 
  geom_boxplot()

# Histogram of all seq_lens
hist(hap$seq_len)

# Genic haplotypes only
genic <- hap[hap$ref_range_id < 35677, ]
ggplot(genic, aes(x = as.factor(gamete_grp_id), y = seq_len)) + 
  geom_boxplot() +
  ylim(0,50000)

# Number of haplotypes per reference range and other basic stats by ref range
temp <- genic %>%
  group_by(ref_range_id) %>%
  summarise(mean = mean(seq_len), min = min(seq_len), max = max(seq_len), n = n())

# Load in count data
consen_count_0.01 <- read.table("Consensus_counts_0.01.txt", header = FALSE, sep = "|", col.names = c("ref_range", "haplotype_count"))
consen_count_0.005 <- read.table("Consensus_counts_0.005.txt", header = FALSE, sep = "|", col.names = c("ref_range", "haplotype_count"))
consen_count_0.001 <- read.table("Consensus_counts_0.001.txt", header = FALSE, sep = ",", col.names = c("ref_range", "haplotype_count"))
consen_count_0.0005 <- read.table("Consensus_counts_0.0005.txt", header = FALSE, sep = "|", col.names = c("ref_range", "haplotype_count"))
consen_count_0.0001 <- read.table("Consensus_counts_0.0001.txt", header = FALSE, sep = ",", col.names = c("ref_range", "haplotype_count"))

consen_count_0.01$annotation <- c(rep("genic", (70754/2)-1), rep("intergenic", 70754/2))
consen_count_0.005$annotation <- c(rep("genic", (70754/2)-1), rep("intergenic", 70754/2))
consen_count_0.001$annotation <-  c(rep("genic", (70754/2)-1), rep("intergenic", 70754/2))
consen_count_0.0005$annotation <- c(rep("genic", (70754/2)-1), rep("intergenic", 70754/2))
consen_count_0.0001$annotation <- c(rep("genic", (70754/2)-1), rep("intergenic", 70754/2))

consen_count_0.01$method <- rep("maxdiv_0.01", nrow(consen_count_0.01))
consen_count_0.005$method <- rep("maxdiv_0.005", nrow(consen_count_0.005))
consen_count_0.001$method <- rep("maxdiv_0.001", nrow(consen_count_0.001))
consen_count_0.0005$method <- rep("maxdiv_0.0005", nrow(consen_count_0.0005))
consen_count_0.0001$method <- rep("maxdiv_0.0001", nrow(consen_count_0.0001))

# Consen count
all_consen_count <- rbind(consen_count_0.0001, consen_count_0.0005, consen_count_0.001, consen_count_0.005, consen_count_0.01)

# Plot consen count by genic and intergenic
ggplot(all_consen_count, aes(x = method, y = haplotype_count, fill = annotation)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_blank(),
        text = element_text(size=20)) +
  labs(x="level of maxDiv collapsing", y = "Number of haplotypes per ref-range") +
  facet_wrap(~annotation, nrow = 2)+ guides(fill=FALSE)

# Plot consen count for both genic and intergenic
ggplot(all_consen_count, aes(x = method, y = haplotype_count)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_blank(),
        text = element_text(size=20)) +
  labs(x="level of maxDiv collapsing", y = "Number of haplotypes per ref-range")

# Ref ranges 
ref <- read.table("reference_ranges.txt", header = F, sep = ",", col.names = c("ref_range", "chrom", "start", "end"))
temp <- merge(x = ref, y= consen_count_0.001, by = "ref_range")
temp <- temp[temp$annotation == "genic",]
# temp$range_len <- temp$end-temp$start
# temp <- temp[temp$chrom==1,]
ggplot(temp, aes(x = start, y = haplotype_count)) + 
  geom_point(aes(colour = factor(annotation))) +
  theme(legend.title = element_blank(),
        text = element_text(size=15)) +
  labs(x="Start Position of Reference Range", y = "Number of haplotypes per ref-range") +
  facet_wrap(~chrom, nrow = 10) +
  guides(fill=FALSE)

# Investigate blank spot on chr5, dominated by large intergenic haplotypes
temp <- temp[which(temp$start > 199000000 & temp$start < 203500000 & temp$chrom == 5), ]


# Subset and investigate regions
temp <- merge(x = ref, y= consen_count_0.001, by = "ref_range")
# temp <- temp[temp$annotation == "genic",]
# temp <- temp[which(temp$start > 124680523 & temp$start < 136656761 & temp$chrom == 8), ]
temp <- temp[which(temp$start > 135576768 & temp$start < 135656761 & temp$chrom == 8), ]

#CCN8, 135653000-135657000
# Rap2.7:135653893 - 13565676
# ZCN8: 126680523 - 126678665

# temp$range_len <- temp$end-temp$start
# temp <- temp[temp$chrom==1,]
ggplot(temp, aes(x = start, y = haplotype_count)) + 
  geom_point(aes(colour = factor(annotation))) +
  theme(legend.title = element_blank(),
        text = element_text(size=15)) +
  labs(x="Start Position of Reference Range", 
       y = "Number of haplotypes per ref-range") +
  ggtitle("Rap2.7 and ZCN8 on Chr 8") +
  facet_wrap(~chrom, nrow = 10) +
  guides(fill=FALSE)

#  su1 locus
# 43430142-43438931
temp <- merge(x = ref, y= consen_count_0.001, by = "ref_range")
su1 <- temp[which(temp$start > 41430142 & temp$start < 45438931 & temp$chrom == 4), ]
su1 <- temp[which(temp$start > 43330142 & temp$start < 43538931 & temp$chrom == 4), ]

ggplot(su1, aes(x = start, y = haplotype_count)) + 
  geom_point(aes(colour = factor(annotation))) +
  theme(legend.title = element_blank(),
        text = element_text(size=15)) +
  labs(x="Start Position of Reference Range", 
       y = "Number of haplotypes per ref-range") +
  ggtitle("su1 - sugary1")

haps <- read.table("ConsensusHapData_001.txt", header = T)
su_haps <- haps[which(haps$chr==4 & haps$start > 43330142 & haps$start < 43538931),]


# -----------------------
# Load in haplotype info
# -----------------------

consen_hap_0.01 <- read.table("Consensus_haplotypes_0.01.txt", header = FALSE, sep = ",", col.names = c("haplotype_id", "gamete_grp_id", "ref_range_id", "seq_len"))
consen_hap_0.005 <- read.table("Consensus_haplotypes_0.005.txt", header = FALSE, sep = ",", col.names = c("haplotype_id", "gamete_grp_id", "ref_range_id", "seq_len"))
consen_hap_0.001 <- read.table("Consensus_haplotypes_0.001.txt", header = FALSE, sep = ",", col.names = c("haplotype_id", "gamete_grp_id", "ref_range_id", "seq_len"))
consen_hap_0.0005 <- read.table("Consensus_haplotypes_0.0005.txt", header = FALSE, sep = ",", col.names = c("haplotype_id", "gamete_grp_id", "ref_range_id", "seq_len"))
consen_hap_0.0001 <- read.table("Consensus_haplotypes_0.0001.txt", header = FALSE, sep = ",", col.names = c("haplotype_id", "gamete_grp_id", "ref_range_id", "seq_len"))

consen_hap_0.01$annotation <- c(rep("genic", (70754/2)-1), rep("intergenic", 70754/2))
consen_hap_0.005$annotation <- c(rep("genic", (70754/2)-1), rep("intergenic", 70754/2))
consen_hap_0.001$annotation <-  c(rep("genic", (70754/2)-1), rep("intergenic", 70754/2))
consen_hap_0.0005$annotation <- c(rep("genic", (70754/2)-1), rep("intergenic", 70754/2))
consen_hap_0.0001$annotation <- c(rep("genic", (70754/2)-1), rep("intergenic", 70754/2))

consen_hap_0.01$method <- rep("maxdiv_0.01", nrow(consen_hap_0.01))
consen_hap_0.005$method <- rep("maxdiv_0.005", nrow(consen_hap_0.005))
consen_hap_0.001$method <- rep("maxdiv_0.001", nrow(consen_hap_0.001))
consen_hap_0.0005$method <- rep("maxdiv_0.0005", nrow(consen_hap_0.0005))
consen_hap_0.0001$method <- rep("maxdiv_0.0001", nrow(consen_hap_0.0001))

# Filter and plot
temp <- consen_count_0.001
temp <- temp[temp$haplotype < 35377, ]
temp <- temp[temp$haplotype < 7000, ]
ggplot(temp, aes(x = haplotype, y = count)) + 
  geom_point()


# -----------------------
# Uncollapsed haplotypes
# -----------------------

# No modification needed
hap <- read.table("cintaMap.txt", header = T) # Lynn generated

# Melt data 
library(reshape)
hap$collapse_level <- rep("uncollapsed", nrow(hap))
hap$annotation <- c(rep("genic", 71354/2), rep("intergenic", 71354/2)) # annotation
mdata <- melt(hap, id=c("refRangeID","chr","start","end", "annotation"))

# Count number of missing haplotypes per taxa
temp <- colSums(hap == 0)

# Plot faceted lengths  of haplotypes by taxa
ggplot(mdata, aes(x = variable, y = log10(value), fill = annotation)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_blank(),
        text = element_text(size=20)) +
  labs(x="Taxa", y = "log10(haplotype length)") +
  facet_wrap(~annotation, nrow = 2)+ guides(fill=FALSE)


# --------------------------------------
#  Investigate aligned haplotype lengths
# --------------------------------------

# Read in file Zack sent 
distances <- read.delim("outputDistances_2ndTry.txt")

# Both haplengths similar
par(mfrow = c(2,1))
hist(distances$taxa1HapLength-distances$denominator)
hist(distances$taxa2HapLength-distances$denominator)

summary(distances$taxa1HapLength-distances$denominator)
summary(distances$taxa2HapLength-distances$denominator)

# Hap1 length against length of denominator
par(mfrow = c(1,1))
plot(log(distances$denominator),log(distances$taxa1HapLength))


# Separate
distances$genic <- distances$refId < max(distances$refId)/2

# Split density bins
ggplot(distances, aes(log10(denominator), log10(taxa1HapLength), fill = log10(..count..)))+
  geom_hex() +
  facet_wrap(~genic) +
  theme(text = element_text(size=20)) +
  labs(x = "log10(denominator len)", x = "log10(Taxa 1 haplotype length)")

# Split histograms
ggplot(distances, aes(distances$denominator/distances$taxa1HapLength))+
  geom_histogram() +
  facet_wrap(~genic) +
  theme(text = element_text(size=20)) +
  labs(x = "log10(denominator len)/log10(Taxa 1 haplotype length)")



