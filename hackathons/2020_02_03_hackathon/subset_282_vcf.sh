#!/bin/bash

# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-02-05
#
# Description 
#   - Gather vcf files with 282
#   - Beagle impute all taxa
#	- Prepare SNPs for two pipelines 1) Calculate PCs, 2) Do GWAS
#	- Both) Filter sites based on imputation accuracy
# ---------------------------------------------------------------


# Subset 282 SNPs
# Subset nam founders
for CHROM in {1..10}
do
  echo "Start data subsetting on chromosome ${CHROM}"
  bcftools index -c /workdir/mbb262/AGPv4/hmp321_282_agpv4_merged_chr${CHROM}.imputed.vcf.gz
  bcftools view -S /workdir/mbb262/goodman282/282set_lines.txt hmp321_282_agpv4_merged_chr${CHROM}.imputed.vcf.gz -Oz -o /workdir/mbb262/goodman282/hmp321_282_agpv4_merged_chr${CHROM}_imputed_goodman282.vcf.gz
  echo "End data subsetting on chromosome ${CHROM}"
done
