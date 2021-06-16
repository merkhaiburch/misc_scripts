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

# Get data from cbsu
# 

# Loop through each chromosome, calculate kinship genome wide
for CHROM in {1..10}
do
  echo "Start analyzing chromosome ${CHROM}"
  /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/ld_matrices/debug_calculate_grm_chr${CHROM}_nam.log \
    -Xmx195g \
    -maxThreads 39 \
    -importGuess /workdir/mbb262/ -noDepth \
    -KinshipPlugin \
    -method, Centered_IBS \
    -endPlugin \
    -export /workdir/mbb262/ \
    -exportType SqrMatrix
  echo "End analyzing chromosome ${CHROM}"
done
