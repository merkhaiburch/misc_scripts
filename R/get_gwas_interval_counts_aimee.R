# -----------------------------------------------------------------------------------------------------------
# Select genes for Aimee
# -----------------------------------------------------------------------------------------------------------

# Note ran script 25 from pleiotropy project to make pleio_ratios_melt3
aimee <- pleio_ratios_melt3 %>% 
  filter(pop %in% c("Goodman Expression", "NAM Field") & range_type == "genic") %>% 
  select("rr_id", "pop", "raw_count")
aimee <- merge(aimee, closest_gene, by = "rr_id") %>% 
  filter(distance == 0) %>% 
  select("pop", "raw_count", "gene_chrom", "gene_start", "gene_end", "gene")
aimee <- merge(aimee, unique(v5), by.x = "gene", by.y = "v4_gene")
dim(aimee)
write.csv(aimee,"~/Downloads/aimee_counts.csv", row.names = FALSE, quote = FALSE)

