#!/usr/bin/env bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-11-15 
#
# Description 
#   - Generate kinship matrices by chromosome for Baoxing
# ---------------------------------------------------------------



# Using NAM imputed Beagle data from Guillaume's directory
# /data1/users/gr226/NAM/Beagle_imputed/imputed
# scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/gr226/NAM/Beagle_imputed/imputed/*.recode.vcf.gz /workdir/mbb262/nam

# Generate a small kinship matrix for testing purposes
for CHROM in {1..10}
do
  echo "Start calculating kinship on chromosome ${CHROM}"
  /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/nam/small_debug_calculate_kinship_nam_chr_${CHROM}.log \
    -Xmx250g \
    -maxThreads 38 \
    -importGuess /workdir/mbb262/nam/AGPv4_parents_chr${CHROM}.recode.vcf.gz -noDepth \
    -FilterSiteBuilderPlugin \
    -siteMinAlleleFreq 0.01 \
    -siteMaxAlleleFreq 0.99 \
    -removeSitesWithIndels true \
    -endPlugin \
    -filterAlign \
    -subsetSites 10000 \
    -step \
    -KinshipPlugin \
    -method Centered_IBS \
    -endPlugin \
    -export /workdir/mbb262/nam/small_kinship_nam_chr_${CHROM}.txt
  echo "Finished calculating kinship on chromosome ${CHROM}"
done


# transfer back to blfs1
# scp /workdir/mbb262/nam/small_kinship_nam_chr_* mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/genotypes/nam/kinship


# Generate the large kinship matrix for baoxing
for CHROM in {1..10}
do
  echo "Start calculating kinship on chromosome ${CHROM}"
  /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/nam/debug_calculate_kinship_nam_chr_${CHROM}.log \
    -Xmx250g \
    -maxThreads 38 \
    -importGuess /workdir/mbb262/nam/AGPv4_parents_chr${CHROM}.recode.vcf.gz -noDepth \
    -FilterSiteBuilderPlugin \
    -siteMinAlleleFreq 0.01 \
    -siteMaxAlleleFreq 0.99 \
    -removeSitesWithIndels true \
    -endPlugin \
    -KinshipPlugin \
    -method Centered_IBS \
    -endPlugin \
    -export /workdir/mbb262/nam/kinship_nam_chr_${CHROM}.txt
  echo "Finished calculating kinship on chromosome ${CHROM}"
done

# transfer back to blfs1
# scp /workdir/mbb262/nam/kinship_nam_chr_* mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/genotypes/nam/kinship
# scp /workdir/mbb262/nam/debug_calculate_kinship_nam_chr_* mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/genotypes/nam/kinship



