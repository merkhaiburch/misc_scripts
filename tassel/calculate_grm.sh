#!/bin/bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-06-09 
#
# Description 
#   - Calculate a genomic relationship matrix for NAM and 282
#   - unimputed SNPs using TASSEL
# ---------------------------------------------------------------

# Get data from cbsu
# scp mbb262@cbsublfs1.biohpc.cornell.edu: /workdir/mbb262


# Loop through each chromosome, calculate LD
for CHROM in {1..10}
do
  echo "Start analyzing chromosome ${CHROM}"
  /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/ld_matrices/debug_calculateLD_chr${CHROM}_nam.log \
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
