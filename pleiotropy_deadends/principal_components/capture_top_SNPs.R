# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-03-09 
#
# Description 
#   - Select top 5% SNPs within global PC windows
#   - Use veriable selection methods to prune for LD, give top SNP
# ---------------------------------------------------------------

# Count number of snps -> 1, 10, 2:9

# IBM Only --> 25M SNPs, n=267
ibm_beagle_imputed <- c(3872582,1849442,2831803,2962916,3259967,2638123,2019003,2188576,2054188,1878419)

#  NAM only --> 25M SNPs, same as IBMs, n=4963
nam_only_imputed <- c(3872582,1849442,2831803,2962916,3259967,2638123,2019003,2188576,2054188,1878419)

# NAM and IBMs, n = 5230, same as above
nam_ibm <- c(3872582,1849442,2831803,2962916,3259967,2638123,2019003,2188576,2054188,1878419)

# After filtering, 1:10, 13M
filtering <- c(2018830, 1431015,1544476,1759736,1336527,1046335,1109426,1024425,903770,984472)

# Model 3 results
m3 <- c(4318192,2066977,2975104,3250270,2942659)

# Load in packages
library(dplyr)

# Load in model 3 DTS results
dts <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2020_02_24_newHapmapSNPs/dts_model3.csv")

# Load in locations of windows
wind <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2020_02_24_newHapmapSNPs/window_file_181genes.csv")

# Load in Z matrix (nxm matrix of taxa by covariates)
Z <- read.table("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/genos/hapmap_v321_snps/ames2nam_3gPCs_allNAM_geneWindow_R.txt", header = TRUE)

# Load in y matrix (nxy matrix of taxa by phenotypes)
y <- read.table("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/genos/hapmap_v321_snps/wallace_phenos_2020_02_24 copy.txt", header = TRUE)
y <- data.frame(y$Trait, y$Days_To_Silk_BLUP_Sum0607_Buckler2009)

# Merge datasets
yZ <- merge(x = y, y = Z, by.x = "y.Trait", by.y = "Trait")


# Load in subsetted X matrix (nxp matrix of taxa by SNPs), substted to only top 5% SNPs
# x.vcf <- "~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2020_02_24_newHapmapSNPs/"
# x_vcf <- snpgdsVCF2GDS(x.vcf, "x.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
# snpgdsSummary("x.gds")
# genofile_x <- snpgdsOpen(x_vcf)

# For loop to capture top 1% of SNPs within a window and add them all to the same dataframe
capture_top_1_snps_window <- function(window_file_df, results){
  top_snp_list <- data.frame()
  for (window in seq(1:nrow(window_file_df))){
    
    # Show me where the loop is
    # print(window_file_df$window[window])
    
    # Subset snps by window
    temp <- results %>% 
      filter(Chr == window_file_df$chrom[window], 
             Pos >= window_file_df$start[window], 
             Pos <= window_file_df$stop[window])
    
    # Calculate a significance threshold for that window
    thresh <- quantile(temp$p, probs = 0.001)
    print(thresh)
    # Subset SNPs in window by p-value
    temp <- temp %>% filter(p <= thresh)
    
    # Add annotation and append to the main dataframe
    temp$info <- rep(window_file_df$window[window], nrow(temp))
    top_snp_list <- rbind(top_snp_list, temp)
  }
  return(top_snp_list)
}


# Get chromosome position file
chr_pos <- top_snp_list[,c(4:5)]
colnames(chr_pos) <- c("Chromosome", "Position")

for (CHROM in seq(1:10)){
  chr <- chr_pos %>% filter(Chromosome == CHROM)
  name <- paste0("chrpos_chr_", CHROM, "_dts.txt")
  write.table(chr, name, 
              row.names = F, col.names = F, quote = F, sep = "\t")
}


# Feed lists of SNPs into bglr?






# Test varbvs
library(varbvs)
# LINEAR REGRESSION EXAMPLE
# -------------------------
# Data are 200 uncorrelated ("unlinked") single nucleotide polymorphisms
# (SNPs) with simulated genotypes, in which the first 20 of them have an
# effect on the outcome. Also generate data for 3 covariates.
maf <- 0.05 + 0.45*runif(200)
X <- (runif(400*200) < maf) + (runif(400*200) < maf)
X <- matrix(as.double(X),400,200,byrow = TRUE)
Z <- randn(400,3)
# Generate the ground-truth regression coefficients for the variables
# (X) and additional 3 covariates (Z). Adjust the QTL effects so that
# the variables (SNPs) explain 50 percent of the variance in the
# outcome.
u <- c(-1,2,1)
beta <- c(rnorm(20),rep(0,180))
beta <- 1/sd(c(X %*% beta)) * beta
# Generate the quantitative trait measurements.
y <- c(-2 + Z %*% u + X %*% beta + rnorm(400))
# Fit the variable selection model.
fit <- varbvs(X,Z,y,logodds = seq(-3,-1,0.1))
print(summary(fit))
# Compute the posterior mean estimate of hyperparameter sa.
sa <- with(fit,sum(sa * w))
# Compare estimated outcomes against observed outcomes.
y.fit <- predict(fit,X,Z)
print(cor(y,y.fit))




# Box 3a: Fitting a model to markers and non-genetic effects in BGLR

#1# Loading and preparing the input data
library(BGLR); data(mice);
Y<-mice.pheno; X<-mice.X; A=mice.A;
y<-Y$Obesity.BMI; y<-(y-mean(y))/sd(y)
#2# Setting the linear predictor
ETA<-list( list(~factor(GENDER)+factor(Litter),
                data=Y,model='FIXED'),
           list(~factor(cage),data=Y, model='BRR'),
           list(X=X, model='BL')
)
#3# Fitting the model
fm<-BGLR(y=y,ETA=ETA, nIter=12000, burnIn=2000)

#1# Estimated Marker Effects & posterior SDs
bHat<- fm$ETA[[3]]$b
SD.bHat<- fm$ETA[[3]]$SD.b
plot(bHat^2, ylab='Estimated Squared-Marker Effect',
     type='o',cex=.5,col=4,main='Marker Effects')

#2# Predictions
# Total prediction
yHat<-fm$yHat
tmp<-range(c(y,yHat))
plot(yHat~y,xlab='Observed',ylab='Predicted',col=2,
     xlim=tmp,ylim=tmp); abline(a=0,b=1,col=4,lwd=2)
# Just the genomic part
gHat<-X%*%fm$ETA[[3]]$b
plot(gHat~y,xlab='Phenotype',
     ylab='Predicted Genomic Value',col=2,
     xlim=tmp,ylim=tmp); abline(a=0,b=1,col=4,lwd=2)
#3# Godness of fit and related statistics
fm$fit
fm$varE # compare to var(y)
#4# Trace plots
list.files()
# Residual variance
varE<-scan('varE.dat')
plot(varE,type='o',col=2,cex=.5,ylab=expression(var[e]));
abline(h=fm$varE,col=4,lwd=2);
abline(v=fm$burnIn/fm$thin,col=4)
# lambda (regularization parameter of the Bayesian Lasso)
lambda<-scan('ETA_3_lambda.dat')
plot(lambda,type='o',col=2,cex=.5,ylab=expression(lambda));
abline(h=fm$ETA[[3]]$lambda,col=4,lwd=2);
abline(v=fm$burnIn/fm$thin,col=4)


# Prior gene list
prior_genes <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/verification_loci/verified_genes_traits.csv",
                        header = T)
prior_genes_10 <- prior_genes %>% filter(seqid == 10, trait == "flowering_time")

# Model 3 just chromosome 10
chr10 <- dts %>% filter(Chr == 10)
source('~/git_projects/haplotype_v_snp_gwas/src/R/random_scripts/manhattan_plot.R')

sig_threshold_1 = -log10(quantile(chr10$p, probs = 0.05))
sig_threshold_2 = -log10(quantile(chr10$p, probs = 0.01))
ylim = 50
title ="Model 3: y ~ SNP + 3gPCs + 660 wPCs for DTS"
# Three windows where four flowering time genes are
temp <- wind[c(211,212,218),]
temp$y1 <- c(0,5,0)
temp$y2 <- c(32,45,45)

ggplot(chr10, aes(x = Pos, y = -log10(p), size = -log10(p))) +
  geom_point(size = 1) +
  geom_hline(yintercept = sig_threshold_1, color = "grey40", linetype = "dashed") +
  geom_hline(yintercept = sig_threshold_2, color = "black", linetype = "dashed") +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "-log10(p)", title = title) + 
  # ggtitle(title) +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 9, vjust = 0.5),
    axis.text.y = element_text(size = 9, vjust = 0.5)) +
  annotate("text", x = 77700138, y = 42, label = "ZmCCA1") +
  geom_vline(xintercept=77700138) +
  
  annotate("text", x = 77722777, y = 40, label = "ZmLHY1") +
  geom_vline(xintercept=77722777) +
  
  annotate("text", x =94430850 , y = 42, label = "ZmCCT") +
  geom_vline(xintercept=94430850) +
  
  annotate("text", x = 141561862, y = 42, label = "ZFL1") +
  geom_vline(xintercept=141561862)+
  
  geom_rect(inherit.aes = F, data=temp, mapping=aes(xmin=start, xmax=stop, ymin=y1, ymax=y2, fill=T), color="black", alpha=0.3)
  

# Model 3 just chromosome 8
chr8 <- dts %>% filter(Chr == 8)
chr8_priors <- prior_genes %>% filter(seqid == 8, trait == "flowering_time")
chr8_priors <- chr8_priors[,4:5]

sig_threshold_1 = -log10(quantile(dts$p, probs = 0.05))
sig_threshold_2 = -log10(quantile(dts$p, probs = 0.01))
ylim = 19
title ="Model 3: y ~ SNP + 3gPCs + 660 wPCs for DTS"
# Three windows where four flowering time genes are
temp <- wind %>% filter(chrom == 8)
temp$y1 <- rep(0, nrow(temp))
temp$y2 <- c(15,18,15,18,15,18,15,18,15,18,15,18,15,18,15,18,15,18,15,18)

ggplot(chr8, aes(x = Pos, y = -log10(p), size = -log10(p))) +
  geom_point(size = 1) +
  geom_hline(yintercept = sig_threshold_1, color = "grey40", linetype = "dashed") +
  geom_hline(yintercept = sig_threshold_2, color = "black", linetype = "dashed") +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "-log10(p)", title = title) + 
  # ggtitle(title) +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 9, vjust = 0.5),
    axis.text.y = element_text(size = 9, vjust = 0.5)) +
  annotate("text", x = 21787238, y = 16, label = "GIGZ1A") +
  geom_vline(xintercept=21787238) +
  
  annotate("text", x = 126880531, y = 16, label = "ZCN8") +
  geom_vline(xintercept=126880531) +
  
  annotate("text", x =136009216 , y = 16.5, label = "ZmRAP2.7") +
  geom_vline(xintercept=136009216) +
  
  annotate("text", x = 163614240, y = 16, label = "ZmHy2") +
  geom_vline(xintercept=163614240)+
  
  geom_rect(inherit.aes = F, data=temp, mapping=aes(xmin=start, xmax=stop, ymin=y1, ymax=y2, fill=T), color="black", alpha=0.3)


