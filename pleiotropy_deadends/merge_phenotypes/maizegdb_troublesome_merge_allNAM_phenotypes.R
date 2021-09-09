# ---------------------------------------
# Found data on MaizeGDB
# https://www.maizegdb.org/traits_ibm_nam
# ---------------------------------------

# Read in data
maize <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/maize_gdb_traits/Trait_Values_IBM_NAM_updated_4-11-16.csv", header = TRUE)

# There are 14 unique studies, of which I did not find data for 5
# Subset out each study, format and merge

# -----------
# Benke 2014
# -----------

benke2014 <- maize %>% filter(Reference == "Benke, A et al. 2014. BMC Plant Biology 14:doi: 10.1186/1471-2229-14-12") %>% select(-Reference, -Stock, -Environment, -Means) %>% droplevels()
benke2014$Trait <- gsub("  Benke 2014", "", benke2014$Trait)
benke2014$Trait <- gsub(", Benke 2014", "", benke2014$Trait)
benke2014$Trait <- gsub(" Benke 2014", "", benke2014$Trait)
benke2014$Trait <- gsub(",", "", benke2014$Trait)
benke2014$Condition[benke2014$Condition == ""] <- NA
levels <- levels(benke2014$Condition)
levels[length(levels) + 1] <- "no_iron_deficiency"
benke2014$Condition <- factor(benke2014$Condition, levels = levels)
benke2014$Condition[is.na(benke2014$Condition)] <- "no_iron_deficiency"
benke2014$Trait <- paste(benke2014$Trait, benke2014$Units, benke2014$Condition,sep = "_")
benke2014$Trait <- gsub(" ", "_", benke2014$Trait)
benke2014 <- benke2014 %>% select(-Plant.Ontology, -Units, -Condition)

# Unmelt data (long to wide format)
benke2014 <- tidyr::pivot_wider(benke2014, names_from = Trait, values_from = Value)

# Save dataset, also gets rid of weird some carryover formatting
write.csv(benke2014, "~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/benke_2014_bmcplantbiology_wholeplant_ironhomeostasis/benke_phenotypes.csv",row.names = FALSE)
benke_2014 <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/benke_2014_bmcplantbiology_wholeplant_ironhomeostasis/benke_phenotypes.csv", header=TRUE)

# Merge with main, weird bugs: 'unknown or uninitialised column' warning
nam_blups <- combine_me_Geno(dataFrame = benke_2014, nameAppend = "raw;Benke_2014", sepBy = "_",
                             byWhat = "Stock.Synonym_raw;Benke_2014")

# ------------
# Burton 2014
# ------------

burton2014 <- maize %>% filter(Reference == "Burton, AL, et al. 2014. Theor Appl Genet. 127:2293-2311")%>% droplevels()

# Clean up trait rows
burton2014$Trait <- gsub("  Burton 2014", "", burton2014$Trait)
burton2014$Trait <- gsub(", Burton 2014", "", burton2014$Trait)
burton2014$Trait <- gsub(" Burton 2014", "", burton2014$Trait)
burton2014$Trait <- gsub(" Burton2014", "", burton2014$Trait)
burton2014$Trait <- gsub(",", "", burton2014$Trait)
burton2014$Trait <- gsub(" ", "_", burton2014$Trait)
burton2014$Trait <- gsub("root__", "root_", burton2014$Trait)
burton2014$Trait <- paste(burton2014$Trait, burton2014$Units, sep = "_")

burton2014 <- burton2014 %>% select(-Plant.Ontology, -Reference, -Stock, -Units, -Means,
                                    -Condition, -Environment)

# Unmelt data (long to wide format)
burton2014 <- tidyr::pivot_wider(burton2014, names_from = Trait, values_from = Value)

# Save dataset, also gets rid of weird some carryover formatting
write.csv(burton2014, "~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/burton_2014_tag_root_rootarchitecture/burton_phenotypes.csv",row.names = FALSE)
burton_2014 <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/burton_2014_tag_root_rootarchitecture/burton_phenotypes.csv", header=TRUE)

# Merge with main, weird bugs: 'unknown or uninitialised column' warning
nam_blups <- combine_me_Geno(dataFrame = burton_2014, nameAppend = "mean;Burton_2014", sepBy = "_",
                             byWhat = "Stock.Synonym_raw;Burton_2014")


# ------------
# Burton 2015
# ------------

burton2015 <- maize %>% filter(Reference == "Burton, AL, et al. 2015. Theor Appl Genet. 128:93-106")%>% droplevels()

burton2015$Trait <- gsub("  Burton 2015", "", burton2015$Trait)
burton2015$Trait <- gsub(", Burton 2015", "", burton2015$Trait)
burton2015$Trait <- gsub(" Burton 2015", "", burton2015$Trait)
burton2015$Trait <- gsub(" Burton2015", "", burton2015$Trait)
burton2015$Trait <- gsub(",", "", burton2015$Trait)
burton2015$Trait <- paste(burton2015$Trait, burton2015$Units, sep = "_")
burton2015$Trait <- gsub(" ", "_", burton2015$Trait)

burton2015 <- burton2015 %>% select(-Plant.Ontology, -Reference, -Stock, -Units, -Means,
                                    -Condition, -Environment)

# Unmelt data (long to wide format)
burton2015 <- tidyr::pivot_wider(burton2015, names_from = Trait, values_from = Value)

# Save dataset, also gets rid of weird some carryover formatting
write.csv(burton2015, "~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/burton_2015_tag_root_rootanatomical/burton_phenotypes.csv",row.names = FALSE)
burton_2015 <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/burton_2015_tag_root_rootanatomical/burton_phenotypes.csv", header=TRUE)

# Merge with main, weird bugs: 'unknown or uninitialised column' warning
nam_blups <- combine_me_Geno(dataFrame = burton_2015, nameAppend = "mean;Burton_2015", sepBy = "_",
                             byWhat = "Stock.Synonym_mean;Burton_2015")
