# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-05-04 
# Updated... 2022-05-04

# Description 
# See if there is an enrichment in GWAS peaks in ID1 binding sites
# and also around ID1 itself
# ---------------------------------------------------------------

# Load in helpful packages
library(dplyr)
library(data.table)
library(ggplot2)
library(GenomicRanges)
library(qqman)

# Load in ID1 binding sites in v4
peaks <- data.table::fread("~/Downloads/ID1_peaksWithSummit.txt")
peaks_sub <- peaks %>% select("Chr", "start", "end", "width")

# Load in flowering time results from buckler 2009 NAM --> NEED TO REMAP AND KEEP P_VALUES
filenames <- list.files("~/Box Sync/Cornell_PhD/labProjects/publications/ID1/data/",
                        pattern = "*filtered.csv",
                        full.names = TRUE)
buckler2009 <- rbindlist(lapply(filenames, data.table::fread))

# Load in flowering time results from buckler 2009 282
filenames <- list.files("/Volumes/merrittData1/pleiotropy/results/goodman_all/",
                        pattern = "*Buckler_2009.txt",
                        full.names = TRUE)
buckler2009 <- rbindlist(lapply(filenames, data.table::fread))

# Subset into dta, dts, and asi results
dta <- buckler2009 %>% filter(Trait == "DTA_BLUP_Buckler_2009_goodman")
dts <- buckler2009 %>% filter(Trait == "DTS_BLUP_Buckler_2009_goodman")
asi <- buckler2009 %>% filter(Trait == "ASI_BLUP_BLUP_Buckler_2009_goodman")

# id1 coordinates: (Chr1: 243199905..243206365), 243,199,905..243,206,365
# id1: Zm00001d032922

# ---------------------------------------------------------------
# Are there any overlaps with id1?
# ---------------------------------------------------------------

# Overlaps of traits within ID1?
temp <- buckler2009 %>% filter(Pos >= 243199905 & Pos <= 243206365 & Chr == 1) # no overlaps

# Overlaps within 5kb of ID1?
temp <- buckler2009 %>% filter(Pos >= 243199905-5000 & Pos <= 243206365+5000 & Chr == 1) # no overlaps

# Overlaps within 10kb of ID1? 
# Need to check if these SNPs are in LD with id1 gene
temp <- buckler2009 %>% filter(Pos >= 243199905-10000 & Pos <= 243206365+10000 & Chr == 1) 

# 4 SNPs, 2 DTA, 2 DTS within 5+ kb downstream
print(temp)


# ---------------------------------------------------------------
# Are there any overlaps with id1 binding sites?
# ---------------------------------------------------------------

# Turn peaks into genomic ranges object
peaks_gr <- GenomicRanges::makeGRangesFromDataFrame(peaks)

# Turn flowering results into object
dta_gr <- GenomicRanges::makeGRangesFromDataFrame(dta,
                                                  seqnames = "Chr", start.field = "Pos", 
                                                  end.field = "Pos", keep.extra.columns=TRUE)
dts_gr <- GenomicRanges::makeGRangesFromDataFrame(dts,
                                                  seqnames = "Chr", start.field = "Pos", 
                                                  end.field = "Pos", keep.extra.columns=TRUE)
asi_gr <- GenomicRanges::makeGRangesFromDataFrame(asi,
                                                  seqnames = "Chr", start.field = "Pos", 
                                                  end.field = "Pos", keep.extra.columns=TRUE)
# Find overlaps between the two
findOverlaps(peaks_gr, dta_gr)
findOverlaps(peaks_gr, dts_gr)
findOverlaps(peaks_gr, asi_gr)


tmp <- subsetByOverlaps(peaks_gr, dta_gr)


# Plot overlaps of TFs and phenotyes ----------------------------------------------------
gwas_df <- dts
gwas_df$p_fdr <- p.adjust(gwas_df$p, method = "BH", n = nrow(gwas_df))

#Get ranges of significant SNPs
significant_snps <- dts_gr
rg_intersect = subsetByOverlaps(significant_snps, peaks_df)
highlight_markers = intersect(gwas_df$Marker, rg_intersect$Marker)

#Read the peak file
all_peaks <- peaks

#Get ranges of peaks for a TF
peaks_df <- peaks_gr

# Make the Manhattan plot on the gwasResults dataset, pretty 
don <- gwas_df %>% 
  group_by(Chr) %>% 
  summarise(chr_len=max(Pos)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(gwas_df, ., by=c("Chr"="Chr")) %>%
  arrange(Chr, Pos) %>%
  mutate(BPcum=Pos+tot) %>%
  mutate(is_highlight=ifelse(Marker %in% highlight_markers, "yes", "no")) %>%
  dplyr::filter(-log10(p)>1.0)

# Make the Manhattan plot on the gwasResults dataset
qqman::manhattan(don, chr="Chr", bp="Pos", snp="Marker", p="p", 
          genomewideline=FALSE, 
          suggestiveline=FALSE,
          highlight=highlight_markers, 
          col= c("grey28","grey61"), 
          main = "DTA GWAS hits overlapping to ID1 binding site peaks")


  
  



