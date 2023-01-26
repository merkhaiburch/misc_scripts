# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-12-28
# Updated... 2022-12-28
#
# Description:
# Test LD regression R package to learn the format of the model
# Manual https://privefl.github.io/bigsnpr/
# ------------------------------------------------------------------------------

# Load packages
install.packages("bigsnpr")
library(bigsnpr)

# Gather variables
ld_score
ld_size
chi2

# Gather fata
bigsnp <- snp_attachExtdata()
G <- bigsnp$genotypes
y <- bigsnp$fam$affection - 1
corr <- snp_cor(G, ind.col = 1:1000)

# Run GWAS model: Column-wise logistic regression
gwas <- big_univLogReg(G, y, ind.col = 1:1000)

# Pull out betas from the model
df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err,
                      n_eff = 4 / (1 / sum(y == 0) + 1 / sum(y == 1)))

# Run the LD regression
snp_ldsc2(corr = corr, df_beta = df_beta)
#>       int        h2 
#> 1.0000000 0.2335429 

# Run LD regression in blocks
snp_ldsc2(corr, df_beta, blocks = 20, intercept = NULL)
