# ---------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-07-19
# Script to merge data together from various
# QTL/gwas studies and appending the author
# name to the trait so they don't get lost
# This script merges all Ames phenotypes together
# ---------------------------------------

ames_phenos <- data.frame()

# -----------
# Romay 2013
# Patient 0
# -----------

romay2013 <- read.table("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/Romay_etal_2013_GenomeBiol/Romay_etal_2013_GenomeBiol_BLUPsphenotypes-130503.txt",
                        header = TRUE, na.strings = "NA")

# Get rid of stuff after semicolon (just isolate line names)
# Looks like sequencing/plate information
romay2013$Complete_name <- gsub(":.*", "", romay2013$Complete_name)

# Add last name
temp <- romay2013[,2:4]
colnames(temp) <- paste(colnames(temp), "BLUP_Romay2013", sep = "_")

# Add to main dataframe
ames_phenos <- cbind(romay2013[,1], temp)

# Edit first column name to something more specific
colnames(ames_phenos)[1] <- "Entry_ID"

# Run consistent name function
ames_phenos$Entry_ID <- gsub("_", "", ames_phenos$Entry_ID)
ames_phenos$Entry_ID <- gsub(" ", "", ames_phenos$Entry_ID)
ames_phenos$Entry_ID <- gsub("-", "", ames_phenos$Entry_ID)
ames_phenos$Entry_ID <- gsub("[.]", "", ames_phenos$Entry_ID)
ames_phenos$Entry_ID <- tolower(ames_phenos$Entry_ID)


# -------------------
# Peiffer et al 2014
# -------------------

# Read data
peiffer2014 <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/Peiffer2014Genetics_blupPhenos20150325.csv",
                        header = TRUE, na.strings = ".")

# Keep only Ames and Assoc panel
peiffer2014 <- peiffer2014[which(peiffer2014$Panel == c("AMES", "ASSO")),]

# Add last name
temp <- peiffer2014[,6:16]
colnames(temp) <- paste(colnames(temp), "BLUP_Peiffer2014", sep = "_")

# Make consistent and simple names for panel, pop, entry, and geno code
# Consistent with Buckler 2006
names <- peiffer2014[,1:5]
colnames(names) <- c("Panel", "pop", "Entry_ID", "entry", "Geno_Code")
names$Entry_ID <- gsub("_", "", names$Entry_ID)
names$Entry_ID <- gsub(" ", "", names$Entry_ID)
names$Entry_ID <- gsub("-", "", names$Entry_ID)
names$Entry_ID <- gsub("[.]", "", names$Entry_ID)
names$Entry_ID <- gsub("\xca", "", names$Entry_ID) # Fix columns with weird encodings
names$Entry_ID <- tolower(names$Entry_ID)

# Add to main data frame
tmp <- cbind(names, temp)
ames_phenos <- merge(x=ames_phenos, y=tmp, 
                     by.x = "Entry_ID", by.y = "Entry_ID", all.x = TRUE)

# Restructure things because it got crazy
ames_phenos <- ames_phenos[,c(1,5:8,2:4,9:19)]


# ----------
# Zila 2014
# ----------

zila2014 <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/Zila2014BMCBiology/FileS1.txt",
                       header = TRUE, sep = ",", na.strings = c(".", "MISSING"))

# Subset data to only DTA, DTS, and rotAVG
zila2014 <- zila2014[,7:10]

# Remove empty rows
zila2014 <- zila2014[rowSums(is.na(zila2014)) != ncol(zila2014),]

# Format names
zila2014$Material <- gsub("_", "", zila2014$Material)
zila2014$Material <- gsub(" ", "", zila2014$Material)
zila2014$Material <- gsub("-", "", zila2014$Material)
zila2014$Material <- gsub("[.]", "", zila2014$Material)
zila2014$Material <- tolower(zila2014$Material)

# Take mean over years/rows
# Take average over rows with same name (average over 3 taxa replicates)
fun1 <- function(x) aggregate(x~Material, zila2014, mean)
dta <- data.frame(apply(zila2014[2], 2, fun1))
dts <- data.frame(apply(zila2014[3], 2, fun1))
rot_AVG <- data.frame(apply(zila2014[4], 2, fun1))

# Merge phenos based on formatted names
tmp <- merge(x=dta, y=dts, by.x = "DTA.Material", by.y = "DTS.Material", all=TRUE)
tmp2 <- merge(x=tmp, y=rot_AVG, by.x = "DTA.Material", by.y = "rot_AVG.Material", all = TRUE)

# Fix column names
names(tmp2) <- gsub(x = names(tmp2), pattern = ".x", replacement = "_mean_raw_Zila2014")
colnames(tmp2)[1] <- "name"

# Merge with main
# Only  1687 /2648 taxa tested 
# Many taxa tested in Zila don't match up with Romay names
ames_phenos <- merge(x=ames_phenos, y=tmp2, 
                     by.x = "Entry_ID", by.y = "name", all.x = TRUE)


# --------------------------------
# Save data in csv and txt format
# --------------------------------

setwd("~/git_projects/haplotype_v_snp_gwas/data")
write.csv(ames_phenos, file="all_Ames_Phenos.csv", row.names = FALSE, quote = FALSE)
write.table(ames_phenos, file = "all_Ames_Phenos.txt", row.names = FALSE, quote = FALSE)




