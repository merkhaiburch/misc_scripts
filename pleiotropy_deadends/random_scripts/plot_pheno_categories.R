# ---------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-06-27
# Script to plot collected pehnotype
# categories in a barchart
# ---------------------------------------

# set directory
setwd("~/Downloads")

# load data
cats <- read.csv("past_mapping_phenotypes - Main_Phenos.csv", header = TRUE)

# Get rid of empty rows
# cats <- cats[c(1:260),]

# Make table
cat_table <- as.data.frame(table(cats$category))

# Change column Name
colnames(cat_table) <- c("trait_category", "frequency")

# Remove weird column at beginning
cat_table <- cat_table[-1,]

# Add an annotation column
cat_table$n <- cat_table$frequency

# Add metabolite column data (dummy data = 100)
# Also remove leaf_metabolite column with only 10 occurances and add here
tmp <- data.frame(as.factor("leaf_metabolites"), 150, 3873+10)
colnames(tmp) <- c("trait_category", "frequency", "n")
cat_table <- rbind(cat_table, tmp)
cat_table <- cat_table[-11,]


# Load packages
library(ggplot2)

# plot
ggplot(cat_table, aes(x = reorder(trait_category, +frequency), y =frequency, fill=trait_category)) + 
  theme_minimal() +
  geom_bar(stat = "identity") +
  theme(legend.position="none", 
        axis.text.x = element_text(angle=45, hjust = 0.7, vjust=1, size = 12),
        axis.text.y = element_text(size = 12)) +
  xlab("Trait Category") +
  ggtitle("Categories of collected NAM, Ames, 282, and IBM  traits") +
  geom_text(aes(label = n), position =  position_dodge(0.9), vjust = -0.5, size=5)

