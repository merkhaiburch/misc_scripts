# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-06-30 
#
# Description 
#   - Assign Plant Trait Ontologies onto phenotypes collected in:
#   - NAM, Goodman Panel, more TBA
# ---------------------------------------------------------------

# Load packages
library(dplyr)
library(tidyr)

# Load in trait file from Google Sheet
df <- read.csv("~/Downloads/past_mapping_phenotypes - main_phenotypes.csv", header = TRUE)

# Melt the data frame by the trait column
df_melt <- df %>% mutate(traits = strsplit(as.character(traits), ",")) %>% unnest(traits)

# Remove NAM rows because we just need Goodman data
df_melt <- df_melt %>% filter(!first_author_name %in% c("Zhou", "Kremling"))

# Export file
write.csv(df_melt, "~/Downloads/past_mapping_phenotypes_revised_melted.csv", row.names = F, quote = F)


# ---------------------
# Add trait ontology
# ---------------------

# Break out the trait to it's more simplified from
df_melt$small_trait <- gsub("_BLUP.*|_BLUE.*|_raw.*", "", df_melt$traits)

# Load in trait ontology


# ----------------------------------------
# One off code to add trait categories 
# to phenotypes for askdb
# ----------------------------------------

# Load packages
library(dplyr)
library(tidyr)


# Load in trait file from Google Sheet
df <- read.csv("~/Downloads/past_mapping_phenotypes - main_phenotypes.csv", header = TRUE)

# Melt the data frame by the trait column
df_melt <- df %>% mutate(traits = strsplit(as.character(traits), ",")) %>% unnest(traits)

# Just get trait and traittype columns
df_melt <- df_melt %>% select(traits, trait_type) %>% filter(trait_type != "metabolite")


# Tack on metabolite and expression names and their trait types (metabolite and expresion)

# Get metabolite names
meta <- read.csv("~/git_projects/haplotype_v_snp_gwas/data/zhou_metabolites_Assoc282_Phenos.csv", header = TRUE)
meta_names <- data.frame(colnames(meta))
meta_names <- meta_names$colnames.meta.[-c(1:11)] %>% data.frame() # remove all annotation that is not trait names
colnames(meta_names)[1] <- "traits" # format column name
meta_names$trait_type <- rep("metabolite", nrow(meta_names))

# Rbind to traits
together <- rbind(df_melt, meta_names)

# Save df and send to Brandon
write.csv(together, file ="~/Downloads/physiological_metabolite_db_traits.csv", quote = F, row.names = F)


# Create same thing for expression datasets --> get summary stat file and get trait names from there

# lead to directory
files <- list.files("Downloads/kremling_expression_v4_pheno_stats/")
setwd("~/Downloads/kremling_expression_v4_pheno_stats/")

# Load in data
expression_dat <- data.frame()
temp_data <- lapply(files, data.table::fread, select = "trait")
expression_dat <- do.call("rbind", temp_data)

# Add in expression name in neighboring column
expression_dat$trait_type <- rep("expression", nrow(expression_dat))
colnames(expression_dat)[1] <- "traits"

# write to file
write.csv(expression_dat, 
          file ="~/Downloads/expression_db_traits.csv", 
          quote = F, row.names = F)


# Format data for google sheet
# Load in trait file from Google Sheet
df <- read.csv("~/Downloads/past_mapping_phenotypes - main_phenotypes.csv", header = TRUE)
df_melt <- df %>% mutate(traits = strsplit(as.character(traits), ",")) %>% unnest(traits)

metabolite_meta <- df_melt %>% filter(trait_type == "metabolite") %>% select(-trait_type, -traits)
tmp <- cbind(metabolite_meta, meta_names)

tmp <- rbind(df_melt, tmp)

# Add in expression data
root <- cbind((expression_dat %>% filter(grepl("GRoot", traits))), df %>% filter(tissue == "germinating root", first_author_name == "Kremling") %>% select(-trait_type, -traits))
shoot <- cbind((expression_dat %>% filter(grepl("GShoot", traits))), df %>% filter(tissue == "germinating shoot", first_author_name == "Kremling") %>% select(-trait_type, -traits))
kern <- cbind((expression_dat %>% filter(grepl("Kern", traits))), df %>% filter(tissue == "kernel", first_author_name == "Kremling") %>% select(-trait_type, -traits))
l3base <- cbind((expression_dat %>% filter(grepl("L3Base", traits))), df %>% filter(tissue == "L3 base", first_author_name == "Kremling") %>% select(-trait_type, -traits))
l3tip <- cbind((expression_dat %>% filter(grepl("L3Tip", traits))), df %>% filter(tissue == "L3 tip", first_author_name == "Kremling") %>% select(-trait_type, -traits))
lman <- cbind((expression_dat %>% filter(grepl("LMAD", traits))), df %>% filter(tissue == "LMAD", first_author_name == "Kremling") %>% select(-trait_type, -traits))
lmad <- cbind((expression_dat %>% filter(grepl("LMAN", traits))), df %>% filter(tissue == "LMAN", first_author_name == "Kremling") %>% select(-trait_type, -traits))

# Rbind everything together
temp <- rbind(tmp, root, shoot, kern, l3base, l3tip, lman, lmad)

# Reorder columns
temp <- temp %>% select(first_author_name, year, species, population, tissue, data_type,
                        publication_name, trait_type, traits, publication_doi, data_doi, notes)

write.csv(temp, "~/Downloads/melted_all_traits_trait_type.csv", quote = F, row.names = F)
