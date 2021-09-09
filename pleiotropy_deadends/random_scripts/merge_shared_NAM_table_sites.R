# ----------------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-07-26
#
# Script to merge table sites in NAM
# that  have a site minimum allele freqeuncy
# of 0.35 (as reccommended) on an individual 
# family-basis, final goal to load into TASSEL
# ----------------------------------------------

# Set directory
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/ames_tests/July_18/FilterMAF035/table2")

# Load all files at once
files <- list.files(pattern = "FilterMAF035_PositionList.txt")
my.data <- lapply(files, read.table, header = TRUE)

# Load IBM file
ibm <- read.table("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2020_01_30_ames2nam_newMarkers/ibmsites_maf0.35_mostdata.txt",
                  header = TRUE)

# Combine lists
tmp <- do.call(c, list(my.data, list(ibm)))

# Isolate lists, turn into dataframes, rbind together
holder <- data.frame()
together <- function(whichList){
  asDF <- data.frame(whichList)
  holder <- rbind(asDF, holder)
  return(holder)
}

# Function to iterate through lists and union join them together
for (i in 1:length(tmp)){
  holder <- together(whichList = tmp[i])
}

# Remove site # column that is making this confusing
holder2 <- holder[,-1]

# Only take unique sites
uniqueSites <- holder2[!duplicated(holder2),]

# Add new site number
tmp <- uniqueSites[with(uniqueSites, order(Chromosome, Position)),]
tmp$Site <- seq(1:nrow(tmp))

# Reorganize
tmp <- tmp[,c(6,1:5)]

# Set directory
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2020_01_30_ames2nam_newMarkers")

# Save this dataframe
write.table(tmp, file = "NAMFamilies_withIBM_unionMarkerJoin0.35.txt", quote = FALSE, row.names = FALSE, sep = '\t')

# Read that dataframe
unionJoinMarkers <- read.table("NAMFamilies_withIBM_unionMarkerJoin0.35.txt", header = TRUE)

# Save as chromosome position file
chrPos <- unionJoinMarkers[,c(3,4)]
write.table(chrPos, file = "NAMFamilies_withIBM_unionMarkerJoin0.35_ChrPos.txt", quote = FALSE, row.names = FALSE, sep = '\t')


# Load IBM sites


