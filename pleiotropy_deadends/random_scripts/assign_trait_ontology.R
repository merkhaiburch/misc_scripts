# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-11-05
#
# Description 
#   - Filtering ontologies from plantemome that were found using Neo4j
#   - Used owl file from here: https://github.com/Planteome/plant-trait-ontology
# ---------------------------------------------------------------

# ----------------------------
# Code used in Neo4j Desktop
# ----------------------------

# Tutorial from Neosemantics: https://neo4j.com/labs/neosemantics/4.0/introduction/

# Neosemantics downloaded here for Neo4j Desktop
# https://install.graphapp.io/

# Read in file directly from github
# CALL n10s.onto.import.fetch("https://raw.githubusercontent.com/Planteome/plant-trait-ontology/master/to.owl", "RDF/XML")

# See subset of all nodes and relationships
# MATCH (n) RETURN *

# Return last ontology nodes and names
# MATCH p=()-->()
# RETURN DISTINCT LAST(nodes(p)).rdfs__label AS ID, LAST(nodes(p)).name AS NAME
# Saved file as a csv


# --------------
# Setup
# --------------

# Load packages
library(dplyr)
library(tidyr)

# Load in file
to <- read.csv("~/Downloads/to_name.csv")

# Neo4j parsed all end nodes, we just want the nodes cooresponding to trait ontology
to <- to %>% filter(grepl("^TO_*", NAME))

# Save file
write.csv(to, "~/git_projects/haplotype_v_snp_gwas/data/trait_ontology_planteome_2020_11_05.csv", row.names = F)


# ----------------------------------------------------
# Formatting trait files and creating melted metadata
# ----------------------------------------------------

# Load in file from google sheets
traits <- read.csv("~/git_projects/haplotype_v_snp_gwas/data/metadata/past_mapping_phenotypes_2020_11_05.csv")

# Melt data on trait column
traits <- traits %>% mutate(traits = strsplit(as.character(traits), ",")) %>% unnest(traits)

# Save file
write.csv(traits, 
          "~/git_projects/haplotype_v_snp_gwas/data/metadata/past_mapping_phenotypes_ontology_2020_11_05.csv", 
          row.names = FALSE)

# Subset by population
nam_traits <- traits %>% filter(population == "NAM")
goodman_traits <- traits %>% filter(population == "Goodman_Association")

# Assign trait ontology to metadata files for all phenotypes







