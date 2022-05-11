# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-04-20 
# Updated... 2022-04-20

# Description 
# Gather gene names and tissues from Karl's expression data
# gene names are in v4, vst transformation
# ---------------------------------------------------------------

# Load in packages
library(dplyr)
library(data.table)

# Gather things outside of the loop
k_names <- c()

# Find all of Karl's data
k_files <- list.files("/Volumes/merrittData1/pleiotropy/phenotypes/kremling", full.names = TRUE)
k_files <- k_files[-4] # get rid of L3

# Gather just tissue names
k_tissues <- gsub("/Volumes/merrittData1/pleiotropy/phenotypes/kremling/|_kremling_formatted_v4_hapmapids.csv", "", k_files)

# For loop iterates though each dataset
for (i in 1:length(k_files)){
  # Load in file
  df <- data.table::fread(k_files[i], nThread = 3)
  
  # gather names
  tmp <- data.frame(colnames(df[,-1]))
  
  # Add tissue identifier
  tmp$tissue <- rep(k_tissues[i], nrow(tmp))
  
  # Change colnames
  colnames(tmp)[1] <- "v4"
  
  # add names to outside thing
  k_names <- rbind(k_names, tmp)
}

# Save to file
write.csv(k_names, "~/Downloads/kremling_names_arcadio.csv", row.names = FALSE, quote = FALSE)
