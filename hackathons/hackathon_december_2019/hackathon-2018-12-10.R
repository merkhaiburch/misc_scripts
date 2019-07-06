# ----------------------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2018-12-10
# 
# Hackathon December 2018
# Team GWAS
# Group: Terry, Travis, Guillaume, Merritt
#
# Test scripts that test ease of use of SNPrelate
# package (and others by discovery) and time it takes
# to import VCF files and eventually run GWAS
# ----------------------------------------------------

# Start clock for measuring system time
start.time <- Sys.time()

# Import/load required packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPRelate", version = "3.8")
library(gdsfmt)
library(SNPRelate)

# Set working directory
setwd("~/Box Sync/Cornell_PhD/Hackathon")

# Import a single VCF file using package
tasselVCF <- snpgdsVCF2GDS("~/Box Sync/Cornell_PhD/Hackathon/rtassel/data/maize_chr9_10thin40000.recode.vcf",
                           "test1.gds")

# Import multiple VCF files (come back to this later)
# testvcf <- c("/Users/mbb262/Box Sync/Cornell_PhD/Hackathon/arabidopsis_thaliana.vcf.gz",
#                 "/Users/mbb262/Box Sync/Cornell_PhD/Hackathon/brachypodium_distachyon.vcf.gz")
# gdsfilesMB <- c("gcf1.gds", "gcf2.gds")
# tassel2VCF <- snpgdsVCF2GDS(testvcf, gdsfilesMB)

# Open GDS file
tassel_genofile <- snpgdsOpen("test1.gds")

# Turn into genotype matrix
genoIntomatrix <- snpgdsGetGeno(tassel_genofile)


# ------------------------------
# GWAS using different packages
# ------------------------------

# Package
library(rrBLUP)

# Fake data and markers from rrGBLUP example
#random population of 200 lines with 1000 markers
M <- matrix(rep(0,200*1000),1000,200)
for (i in 1:200) {
  M[,i] <- ifelse(runif(1000)<0.5,-1,1)
}
colnames(M) <- 1:200
geno <- data.frame(marker=1:1000,chrom=rep(1,1000),pos=1:1000,M,check.names=FALSE)

QTL <- 100*(1:5) #pick 5 QTL
u <- rep(0,1000) #marker effects
u[QTL] <- 1
g <- as.vector(crossprod(M,u))
h2 <- 0.5
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))
pheno <- data.frame(line=1:200,y=y)

# Run GWAS on this data
scores <- GWAS(pheno,geno,plot=T, n.PC = 1)

# Histogram of phennontypes
hist(pheno$y)

# Close file
# snpgdsClose(tassel_genofile)

# End clock for system time
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


