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
library(phangorn)

# Load in data
skmer100k <- read.table("/workdir/hackathon/skmer_results/jc-dist-mat-pandand16.txt", header = TRUE) %>% as.matrix()
skmer100k <- read.table("~/Downloads/jc-dist-mat-pandand16.txt", header = TRUE) %>% as.matrix()


# Change rownames
rownames(skmer100k) <- skmer100k[,1]
skmer100k <- skmer100k[,-1]

# Format matrix using ape
my_nj <- ape::njs(skmer100k)

# Plot the tree as an “unrooted” tree
plot(my_nj, "phylogram")

# Plot as a rooted tree
my_nj_rooted <- ape::root(my_nj,7:8)
plot(my_nj_rooted, "phylogram")


# --------------------
# Plot 1M read tree
# --------------------

skmer1m <- read.table("~/Downloads/jc-dist-mat-pandand16_1m.txt", header = TRUE) %>% as.matrix()

# Change rownames
rownames(skmer1m) <- skmer1m[,1]
skmer1m <- skmer1m[,-1]

# Format matrix using ape
my_nj_1m <- ape::njs(skmer1m)

# Plot the tree as an “unrooted” tree
plot(my_nj_1m, "phylogram")

# Plot as a rooted tree
my_nj_rooted_1m <- ape::root(my_nj_1m,7:8) # check columns for correct root
plot(my_nj_rooted_1m, "phylogram")



# --------------------
# Plot 1M + mew gemomes read tree
# --------------------

skmer1m_plus_new <- read.table("~/Downloads/panand_skmer_16_genomes_1m_plus_new.txt", header = TRUE) %>% as.matrix()

# Change rownames
rownames(skmer1m_plus_new) <- skmer1m_plus_new[,1]
skmer1m_plus_new <- skmer1m_plus_new[,-1]

# Format matrix using ape
my_nj_1m_plus_new <- ape::njs(skmer1m_plus_new)

# Plot as a rooted tree
my_nj_rooted_1m_plus_new <- ape::root(my_nj_1m_plus_new,7:8) # check columns for correct root
plot(my_nj_rooted_1m_plus_new, "phylogram")


