#!/usr/bin/env bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-11-15 
#
# Description 
#   - Example script to run MLM in Tassel
# ---------------------------------------------------------------

# Using NAM imputed Beagle data from Guillaume's directory
# Example name
# AGPv4_NAM_chr1.imputed.vcf.gz
# scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/gr226/NAM/Beagle_imputed/imputed/*.imputed.vcf.gz /workdir/mbb262/nam


# Subset sites within each vcf file, export
for CHROM in {1..10}
do
  echo "I am on chr ${CHROM}"

  /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/nam/subset_nam_chr${CHROM}.log \
    -Xmx1020g \
    -maxThreads 87 \
    -importGuess /workdir/mbb262/nam/AGPv4_NAM_chr${CHROM}.imputed.vcf.gz -noDepth \
    -FilterSiteBuilderPlugin \
    -siteMinAlleleFreq 0.01 \
    -siteMaxAlleleFreq 0.99 \
    -removeSitesWithIndels true \
    -endPlugin \
    -filterAlign \
    -subsetSites 300000 \
    -step \
    -export /workdir/mbb262/nam/subset_nam_snps_chr${CHROM}.vcf \
    -exportType VCF
  echo "I just finished chr ${CHROM}"
done


# combine subsetted files (see if works better)
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
  -debug /workdir/mbb262/nam/debug_merge_vcf_nam_subset.log \
  -Xmx1020g \
  -maxThreads 87 \
  -fork1 -importGuess /workdir/mbb262/nam/subset_nam_snps_chr1.vcf -noDepth \
  -fork2 -importGuess /workdir/mbb262/nam/subset_nam_snps_chr2.vcf -noDepth \
  -fork3 -importGuess /workdir/mbb262/nam/subset_nam_snps_chr3.vcf -noDepth \
  -fork4 -importGuess /workdir/mbb262/nam/subset_nam_snps_chr4.vcf -noDepth \
  -fork5 -importGuess /workdir/mbb262/nam/subset_nam_snps_chr5.vcf -noDepth \
  -fork6 -importGuess /workdir/mbb262/nam/subset_nam_snps_chr6.vcf -noDepth \
  -fork7 -importGuess /workdir/mbb262/nam/subset_nam_snps_chr7.vcf -noDepth \
  -fork8 -importGuess /workdir/mbb262/nam/subset_nam_snps_chr8.vcf -noDepth \
  -fork9 -importGuess /workdir/mbb262/nam/subset_nam_snps_chr9.vcf -noDepth \
  -fork10 -importGuess /workdir/mbb262/nam/subset_nam_snps_chr10.vcf -noDepth \
  -combine11 -input1 -input2 -input3 -input4 -input5 -input6 -input7 -input8 -input9 -input10 \
  -mergeGenotypeTables \
  -export /workdir/mbb262/nam/subset_nam_combined_all_chroms.vcf \
  -exportType VCF


# zip file
bgzip --threads 87 /workdir/mbb262/nam/subset_nam_combined_all_chroms.vcf


# calculate kinship for all chroms at once on subsetted, combined file
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
  -debug /workdir/mbb262/nam/debug_calculate_kinship_nam_all_chroms_subset.log \
  -Xmx1020g \
  -maxThreads 87 \
  -importGuess /workdir/mbb262/nam/subset_nam_combined_all_chroms.vcf.gz -noDepth \
  -FilterSiteBuilderPlugin \
  -siteMinAlleleFreq 0.01 \
  -siteMaxAlleleFreq 0.99 \
  -removeSitesWithIndels true \
  -endPlugin \
  -KinshipPlugin \
  -method Centered_IBS \
  -endPlugin \
  -export /workdir/mbb262/nam/kinship_nam_all_chroms_subset.txt


# transfer all back to blfs1
# kinship matrix
# scp -r /workdir/mbb262/nam/* mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/genotypes/nam/kinship


# # Example to run MLM with kinship matrix for
# y ~ nam SNP + protein structure as covariates + kinship matrix
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
  -debug /debug_MLM_NAMkinship.log \
  -Xmx250g \
  -maxThreads 38 \
  -importGuess /nam_combined_all_chroms.vcf \
  -r /phenotypes.txt \
  -r /protein_or_other_fixed_effects.txt \
  -k /kinship_nam_all_chroms_subset.txt \
  -mlm \
  -mlmVarCompEst P3D \
  -mlmCompressionLevel None \
  -export /mlm_nam_kinship_results


# Calculate PCs on a vcf file

# Filter vcf file
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
  -debug /debug_pca_vcf_nam_subset.log \
  -Xmx1020g \
  -maxThreads 87 \
  -fork1 -importGuess /subset_nam_combined_all_chroms.vcf.gz -noDepth \
  -FilterSiteBuilderPlugin \
  -siteMinAlleleFreq 0.01 \
  -siteMaxAlleleFreq 0.99 \
  -removeSitesWithIndels true \
  -endPlugin \
  -PrincipalComponentsPlugin \
  -covariance true \
  -endPlugin \
  -export /pc_output
