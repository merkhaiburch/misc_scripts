#!/bin/bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-06-09 
#
# Description 
#   - Calculate a genomic relationship matrix for NAM and 282
#   - Use PHG imputed SNPs for:
#   - 1) genome wide SNPs 2) SNPs only within specific intervals
# ---------------------------------------------------------------

# --------------------
#  NAM SNPs and GRM
# --------------------

# Filter data using tassel using MAF then randomly subset 10M sites genome wide
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/debug_calculate_grm_nam.log \
    -Xmx1000g \
    -maxThreads 86 \
    -importGuess /workdir/mbb262/NAM_phg_snps.vcf -noDepth \
    -FilterSiteBuilderPlugin \
    -siteMinAlleleFreq 0.01 \
    -siteMaxAlleleFreq 0.99 \
    -removeSitesWithIndels true \
    -endPlugin \
    -filterAlign \
    -subsetSites 10000000 \
    -step \
    -export /workdir/mbb262/filtered_NAM_phg_snps.vcf \
    -exportType VCF

bgzip --threads 80 /workdir/mbb262/filtered_NAM_phg_snps.vcf

# calculate kinship genome wide using a subsetted vcf
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/debug_calculate_grm_kinship_nam.log \
    -Xmx1000g \
    -maxThreads 86 \
    -importGuess /workdir/mbb262/filtered_NAM_phg_snps.vcf.gz -noDepth \
    -KinshipPlugin \
    -method Centered_IBS \
    -endPlugin \
    -export /workdir/mbb262/kinship_filtered_NAM_phg_snps.txt



# -----------------------------------------
#  NAM SNPs and GRM --> HARE regions only
# -----------------------------------------

# bed file from HARE genes in v3, filtering SNPs within these intervals
# Filter data using tassel using MAF 
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/nam_hare_downsample_debug.log \
    -Xmx1000g \
    -maxThreads 86 \
    -importGuess /workdir/mbb262/NAM_phg_snps.vcf -noDepth \
    -FilterSiteBuilderPlugin \
    -bedFile ./hare_v5_intervals.bed \
    -endPlugin \
    -export /workdir/mbb262/filtered_hare_nam_phg_snps.vcf \
    -exportType VCF

bgzip --threads 80 /workdir/mbb262/filtered_hare_nam_phg_snps.vcf

# calculate kinship genome wide using a subsetted vcf
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/nam_hare_calc_grm_debug.log \
    -Xmx1000g \
    -maxThreads 86 \
    -importGuess /workdir/mbb262/filtered_hare_nam_phg_snps.vcf.gz -noDepth \
    -KinshipPlugin \
    -method Centered_IBS \
    -endPlugin \
    -export /workdir/mbb262/kinship_filtered_nam_hare_phg_snps.txt


# --------------------
#  282 SNPs filter and calc GRM
# --------------------

# Filter data using tassel using MAF 
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/goodman_downsample_debug.log \
    -Xmx500g \
    -maxThreads 60 \
    -importGuess /workdir/ag2484/goodman_phg_snps.txt.vcf -noDepth \
    -filterAlign \
    -filterAlignMinFreq 0.01 \
    -filterAlignMaxFreq 0.99 \
    -export /workdir/mbb262/filtered_goodman_phg_snps.vcf \
    -exportType VCF

bgzip --threads 60 /workdir/mbb262/filtered_goodman_phg_snps.vcf

# calculate kinship genome wide using a subsetted vcf
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/goodman_calc_grm_debug.log \
    -Xmx500g \
    -maxThreads 60 \
    -importGuess /workdir/mbb262/filtered_goodman_phg_snps.vcf.gz -noDepth \
    -KinshipPlugin \
    -method Centered_IBS \
    -endPlugin \
    -export /workdir/mbb262/kinship_filtered_goodman_phg_snps.txt


# -----------------------------------------
#  282 SNPs, filter, calc GRM --> HARE regions only
# -----------------------------------------

# Filter data using tassel using MAF 
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/goodman_hare_downsample_debug.log \
    -Xmx500g \
    -maxThreads 60 \
    -importGuess /workdir/ag2484/goodman_phg_snps.txt.vcf -noDepth \
    -FilterSiteBuilderPlugin \
    -siteMinAlleleFreq 0.01 \
    -siteMaxAlleleFreq 0.99 \
    -removeSitesWithIndels true \
    -bedFile hare_v5_intervals.bed \
    -endPlugin \
    -export /workdir/mbb262/filtered_hare_goodman_phg_snps.vcf \
    -exportType VCF

bgzip --threads 60 /workdir/mbb262/filtered_hare_goodman_phg_snps.vcf

# calculate kinship genome wide using a subsetted vcf
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/complete_grms_goodman/goodman_hare_calc_grm_debug.log \
    -Xmx500g \
    -maxThreads 60 \
    -importGuess /workdir/mbb262/filtered_hare_goodman_phg_snps.vcf.gz -noDepth \
    -KinshipPlugin \
    -method Centered_IBS \
    -endPlugin \
    -export /workdir/mbb262/complete_grms_goodman/kinship_filtered_goodman_hare_phg_snps.txt




