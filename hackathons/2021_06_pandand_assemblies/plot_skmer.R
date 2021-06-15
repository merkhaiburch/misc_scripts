# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-06-15 
#
# Description 
#   - Plot Skmer distance matrix trees
#   - Example: https://rpubs.com/WalshJake75/674724
# ---------------------------------------------------------------

# Load packages
library(ape)
library(magrittr)

# Load in data
skmer100k <- read.table("/workdir/hackathon/skmer_results/jc-dist-mat-pandand16.txt", header = TRUE) %>% as.matrix()

# Change rownames
rownames(skmer100k) <- skmer100k[,1]
skmer100k <- skmer100k[,-1]

# Format matrix using ape
my_nj <- ape::njs(skmer100k)

# Plot the tree as an “unrooted” tree
plot(my_nj, "phylogram")
