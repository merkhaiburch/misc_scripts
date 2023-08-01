#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-07-28
# Updated... 2023-07-28
#
# Description:
# Subsample imputed NAM hybrid vcfs, merge, calculate MDS for Michelle
# NOTE: the DTS_taxa_PHZ51_hybrid.txt and GY_taxa_PHZ51_hybrid.txt are taxa lists from
#       Michelle, but in R I added on "/PHZ51" to the end of the taxa names to
#       be able to filter Guillaume's hybrid vcf files (i.e. Z)
# ------------------------------------------------------------------------------

# Make a direcotry
mkdir -p /workdir/mbb262/mds/genotypes
mkdir -p /workdir/mbb262/mds/subsampled_genotypes_DTS
mkdir -p /workdir/mbb262/mds/subsampled_genotypes_GY

# Get files from Guillaume's directory
scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/gr226/Ramstein_AmesNAMHybrids_2019/NAM/genotypes/hybrid_imputed/*.gz /workdir/mbb262/mds/genotypes

# Subset vcf files with a taxa list, max AF, then by number of sites for DTS
for FILE in *.gz
do
  echo "I am on file ${FILE}"

  /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/mds/subsampled_genotypes_DTS/debug_chr${FILE}.log \
    -Xmx200g \
    -maxThreads 38 \
    -importGuess ./$FILE \
    -noDepth \
    -FilterTaxaBuilderPlugin \
    -taxaList /workdir/mbb262/mds/DTS_taxa_PHZ51_hybrid.txt \
    -endPlugin \
    -filterAlign \
    -subsetSites 400000 \
    -step \
    -export /workdir/mbb262/mds/subsampled_genotypes_DTS/${FILE}_subsampled_DTS.vcf \
    -exportType VCF
  echo "I just finished chr ${FILE}"
done

# Subset vcf files with a taxa list, max AF, then by number of sites for GY
for FILE in *.gz
do
  echo "I am on file ${FILE}"

  /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/mds/subsampled_genotypes_GY/debug_chr${FILE}.log \
    -Xmx200g \
    -maxThreads 38 \
    -importGuess ./$FILE \
    -noDepth \
    -FilterTaxaBuilderPlugin \
    -taxaList /workdir/mbb262/mds/GY_taxa_PHZ51_hybrid.txt \
    -endPlugin \
    -FilterSiteBuilderPlugin \
    -siteMaxAlleleFreq 1.0 \
    -removeSitesWithIndels true \
    -endPlugin \
    -filterAlign \
    -subsetSites 400000 \
    -step \
    -export /workdir/mbb262/mds/subsampled_genotypes_GY/${FILE}_subsampled_GY.vcf \
    -exportType VCF
  echo "I just finished chr ${FILE}"
done


# combine subsetted files for DTS
VCF_PATH=/workdir/mbb262/mds/subsampled_genotypes_DTS
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
  -debug /workdir/mbb262/mds/subsampled_genotypes_DTS/debug_merge_vcf_nam_subset.log \
  -Xmx200g \
  -maxThreads 39 \
  -fork1 -importGuess $VCF_PATH/AGPv4_NAM_chr1.PHZ51-het.vcf.gz_subsampled_DTS.vcf -noDepth \
  -fork2 -importGuess $VCF_PATH/AGPv4_NAM_chr2.PHZ51-het.vcf.gz_subsampled_DTS.vcf -noDepth \
  -fork3 -importGuess $VCF_PATH/AGPv4_NAM_chr3.PHZ51-het.vcf.gz_subsampled_DTS.vcf -noDepth \
  -fork4 -importGuess $VCF_PATH/AGPv4_NAM_chr4.PHZ51-het.vcf.gz_subsampled_DTS.vcf -noDepth \
  -fork5 -importGuess $VCF_PATH/AGPv4_NAM_chr5.PHZ51-het.vcf.gz_subsampled_DTS.vcf -noDepth \
  -fork6 -importGuess $VCF_PATH/AGPv4_NAM_chr6.PHZ51-het.vcf.gz_subsampled_DTS.vcf -noDepth \
  -fork7 -importGuess $VCF_PATH/AGPv4_NAM_chr7.PHZ51-het.vcf.gz_subsampled_DTS.vcf -noDepth \
  -fork8 -importGuess $VCF_PATH/AGPv4_NAM_chr8.PHZ51-het.vcf.gz_subsampled_DTS.vcf -noDepth \
  -fork9 -importGuess $VCF_PATH/AGPv4_NAM_chr9.PHZ51-het.vcf.gz_subsampled_DTS.vcf -noDepth \
  -fork10 -importGuess $VCF_PATH/AGPv4_NAM_chr10.PHZ51-het.vcf.gz_subsampled_DTS.vcf -noDepth \
  -combine11 -input1 -input2 -input3 -input4 -input5 -input6 -input7 -input8 -input9 -input10 \
  -mergeGenotypeTables \
  -export $VCF_PATH/AGPv4_NAM_chrALL.PHZ51-het_DTS.vcf \
  -exportType VCF
gzip $VCF_PATH/AGPv4_NAM_chrALL.PHZ51-het_DTS.vcf

# combine subsetted files for GY
VCF_PATH=/workdir/mbb262/mds/subsampled_genotypes_GY
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
  -debug /workdir/mbb262/mds/subsampled_genotypes_GY/debug_merge_vcf_nam_subset.log \
  -Xmx200g \
  -maxThreads 39 \
  -fork1 -importGuess $VCF_PATH/AGPv4_NAM_chr1.PHZ51-het.vcf.gz_subsampled_GY.vcf -noDepth \
  -fork2 -importGuess $VCF_PATH/AGPv4_NAM_chr2.PHZ51-het.vcf.gz_subsampled_GY.vcf -noDepth \
  -fork3 -importGuess $VCF_PATH/AGPv4_NAM_chr3.PHZ51-het.vcf.gz_subsampled_GY.vcf -noDepth \
  -fork4 -importGuess $VCF_PATH/AGPv4_NAM_chr4.PHZ51-het.vcf.gz_subsampled_GY.vcf -noDepth \
  -fork5 -importGuess $VCF_PATH/AGPv4_NAM_chr5.PHZ51-het.vcf.gz_subsampled_GY.vcf -noDepth \
  -fork6 -importGuess $VCF_PATH/AGPv4_NAM_chr6.PHZ51-het.vcf.gz_subsampled_GY.vcf -noDepth \
  -fork7 -importGuess $VCF_PATH/AGPv4_NAM_chr7.PHZ51-het.vcf.gz_subsampled_GY.vcf -noDepth \
  -fork8 -importGuess $VCF_PATH/AGPv4_NAM_chr8.PHZ51-het.vcf.gz_subsampled_GY.vcf -noDepth \
  -fork9 -importGuess $VCF_PATH/AGPv4_NAM_chr9.PHZ51-het.vcf.gz_subsampled_GY.vcf -noDepth \
  -fork10 -importGuess $VCF_PATH/AGPv4_NAM_chr10.PHZ51-het.vcf.gz_subsampled_GY.vcf -noDepth \
  -combine11 -input1 -input2 -input3 -input4 -input5 -input6 -input7 -input8 -input9 -input10 \
  -mergeGenotypeTables \
  -export $VCF_PATH/AGPv4_NAM_chrALL.PHZ51-het_GY.vcf \
  -exportType VCF
gzip $VCF_PATH/AGPv4_NAM_chrALL.PHZ51-het_GY.vcf


# --------------------------------------------------------------------------------
# Calculate a distance matrix for each vcf file
# --------------------------------------------------------------------------------

# DTS
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/mds/DTS_debug_distance.log \
    -Xmx200g \
    -maxThreads 38 \
    -importGuess /workdir/mbb262/mds/subsampled_genotypes_DTS/AGPv4_NAM_chrALL.PHZ51-het_DTS.vcf.gz \
    -noDepth \
    -DistanceMatrixPlugin \
    -endPlugin \
    -MultiDimensionalScalingPlugin \
    -endPlugin \
    -export /workdir/mbb262/mds/mds_DTS.txt



# GY
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/mds/GY_debug_distance.log \
    -Xmx200g \
    -maxThreads 38 \
    -importGuess /workdir/mbb262/mds/subsampled_genotypes_GY/AGPv4_NAM_chrALL.PHZ51-het_GY.vcf.gz \
    -noDepth \
    -DistanceMatrixPlugin \
    -endPlugin \
    -MultiDimensionalScalingPlugin \
    -endPlugin \
    -export /workdir/mbb262/mds/mds_GY.txt






