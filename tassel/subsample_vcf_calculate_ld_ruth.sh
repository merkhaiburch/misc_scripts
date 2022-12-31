#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-10-05
# Updated... 2022-10-05
#
# Description:
# Subsample SNPs by chromosome for Ruth, calculate LD for each chromosome
# ------------------------------------------------------------------------------

for CHROM in {1..10}
do
  echo "Start subsample SNPs for chromosome ${CHROM}"
  /tassel-5-standalone/run_pipeline.pl \
    -debug ./debug_subset_sites.log \
    -Xmx100g \
    -maxThreads 30 \
    -importGuess ./your_vcf_file.vcf -noDepth \
    -FilterSiteBuilderPlugin \
    -siteMinAlleleFreq 0.01 \
    -siteMaxAlleleFreq 0.99 \
    -startChr $CHROM \
    -endChr $CHROM \
    -endPlugin \
    -filterAlign \
    -subsetSites 1000000 \
    -step \
    -export ./filtered_snps_chr_${CHROM}.vcf \
    -exportType VCF

    echo "Finished subsample SNPs for chromosome ${CHROM}"
done

# zip all the files
gzip *.vcf


# Calculate LD for each chromosome
for CHROM in {1..10}
do
  echo "Start calculating LD on chromosome ${CHROM}"
  /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug ./nam/debug_chr_${CHROM}.log \
    -Xmx250g \
    -maxThreads 30 \
    -importGuess ./your_vcf_chr${CHROM}.vcf.gz -noDepth \
    -FilterSiteBuilderPlugin \
    -siteMinAlleleFreq 0.01 \
    -siteMaxAlleleFreq 0.99 \
    -endPlugin \
    -ld \
    -ldType All \
    -export ./ld_chr_${CHROM}.txt
  echo "Finished calculating LD on chromosome ${CHROM}"
done


