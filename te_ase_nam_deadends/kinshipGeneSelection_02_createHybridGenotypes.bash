#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-07-28 
# Updated... 2022-07-28

# Description 
# Create hybrid genotypes from inbred names
# ---------------------------------------------------------------

# Get v4 --> v5 uplifted data from blfs1 (did uplifting in uplift_hapmap_v5)
mkdir -p /workdir/mbb262/genotypes/inbred_v5
# scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/genotypes/Hapmap321_Ames_v5_uplifted/Hapmap321_inbreds_v5_uplifted/* /workdir/mbb262/genotypes/inbred_v5

# inbred combinations were created in the subset_hapmap_taxa.R script
# the "te_ase_hybrids_for_tassel.txt" file

# Variables
mkdir -p /workdir/mbb262/genotypes/v5_hybrid_genotypes
VCF_OUT=/workdir/mbb262/genotypes/v5_hybrid_genotypes
VCF_IN=/workdir/mbb262/genotypes/inbred_v5
VCF_OUTMERGED=/workdir/mbb262/genotypes/v5_inbred_hybrid_genotypes


# Sort vcf files
for CHROM in {1..10}
do
  /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug $VCF_IN/debug_sorting_chr${CHROM}.log \
    -Xmx200g \
    -maxThreads 30 \
    -importGuess $VCF_IN/hmp321_282_agpv5_uplifted_merged_chr${CHROM}_imputed.vcf.gz \
    -sortPositions \
    -export $VCF_IN/hmp321_282_agpv5_uplifted_merged_chr${CHROM}_imputed_sorted.vcf.gz \
    -exportType VCF
done


# see if there are a similar number of snps in the sorted and unsorted files
# THERE IS AN EQUAL NUMBER
# for CHROM in {1..10}
# do
#   echo "I am on chr ${CHROM}"
#   grep -v "^#" ./hmp321_282_agpv5_uplifted_merged_chr${CHROM}_imputed.vcf | wc -l
#   grep -v "^#" ./hmp321_282_agpv5_uplifted_merged_chr${CHROM}_imputed_sorted.vcf | wc -l
# done


# Run tassel to create hybrid genotypes for:
# B73 x NAM inbreds, PHZ51 x NAM inbreds, and PHB47 x NAM inbreds
for CHROM in {1..10}
do
  /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug $VCF_OUT/debug_hybrid_chr${CHROM}.log \
    -Xmx200g \
    -maxThreads 30 \
    -importGuess $VCF_IN/hmp321_282_agpv5_uplifted_merged_chr${CHROM}_imputed_sorted.vcf.gz \
    -noDepth \
    -CreateHybridGenotypesPlugin \
    -hybridFile /workdir/mbb262/genotypes/te_ase_hybrids_for_tassel.txt \
    -hybridChar "_" \
    -endPlugin \
    -export $VCF_OUT/hmp321_282_agpv5_hybrid_b73_phz51_phb47_uplifted_merged_chr${CHROM}_imputed_sorted \
    -exportType VCF

  bgzip --threads 30 $VCF_OUT/hmp321_282_agpv5_hybrid_b73_phz51_phb47_uplifted_merged_chr${CHROM}_imputed_sorted
done


# Get all taxa names to add to subsetting file
bcftools query -l hmp321_282_agpv5_hybrid_b73_phz51_phb47_uplifted_merged_chr10_imputed_sorted.vcf.gz


# Merge inbred and hybrid vcf files by chromosome
# remove taxa that aren't needed based on txt file
mkdir -p /workdir/mbb262/genotypes/v5_inbred_hybrid_genotypes
for CHROM in {1..10}
do
  echo "I am on chr ${CHROM}"
  /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug $VCF_OUTMERGED/debug_merge_vcf_nam_subset_chr${CHROM}.log \
    -Xmx200g \
    -maxThreads 35 \
    -fork1 -importGuess $VCF_OUT/hmp321_282_agpv5_hybrid_b73_phz51_phb47_uplifted_merged_chr${CHROM}_imputed_sorted.vcf.gz -noDepth \
    -fork2 -importGuess $VCF_IN/hmp321_282_agpv5_uplifted_merged_chr${CHROM}_imputed_sorted.vcf.gz -noDepth \
    -combine3 -input1 -input2 \
    -MergeGenotypeTablesPlugin \
    -endPlugin \
    -FilterSiteBuilderPlugin \
    -startChr ${CHROM} \
    -endChr ${CHROM} \
    -endPlugin \
    -FilterTaxaBuilderPlugin \
    -taxaList /workdir/mbb262/te_ase_nam_inbreds_hybrids_hapmapids.txt\
    -endPlugin \
    -export $VCF_OUTMERGED/hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr${CHROM}_imputed_sorted \
    -exportType VCF

  bgzip --threads 30 $VCF_OUTMERGED/hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr${CHROM}_imputed_sorted.vcf

  echo "I just finished chr ${CHROM}"
done



# Version 2
for CHROM in {1..10}
do
  echo "I am on chr ${CHROM}"
  /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug $VCF_OUTMERGED/debug_merge_vcf_nam_subset.log \
    -Xmx200g \
    -maxThreads 35 \
    -fork1 -importGuess $VCF_OUT/hmp321_282_agpv5_hybrid_b73_phz51_phb47_uplifted_merged_chr10_imputed_sorted.vcf.gz -noDepth \
    -fork2 -importGuess $VCF_IN/hmp321_282_agpv5_uplifted_merged_chr10_imputed_sorted.vcf.gz -noDepth \
    -combine3 -input1 -input2 \
    -MergeGenotypeTablesPlugin \
    -FilterSiteBuilderPlugin \
    -startChr 10 \
    -endChr 10 \
    -endPlugin \
    -FilterTaxaBuilderPlugin \
    -taxaList /workdir/mbb262/te_ase_nam_inbreds_hybrids_hapmapids.txt\
    -endPlugin \
    -export $VCF_OUTMERGED/hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr10_imputed_sorted \
    -exportType VCF

  bgzip --threads 30 $VCF_OUTMERGED/hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr10_imputed_sorted.vcf

  echo "I just finished chr ${CHROM}"
done





# Subsample combined hybrid and inbred vcf file
for CHROM in {1..10}
do
  echo "I am on chr ${CHROM}"

  /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug $VCF_OUTMERGED/subset_nam_chr${CHROM}.log \
    -Xmx200g \
    -maxThreads 30 \
    -importGuess $VCF_OUTMERGED/hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr${CHROM}_imputed_sorted.vcf.gz -noDepth \
    -FilterSiteBuilderPlugin \
    -siteMinAlleleFreq 0.01 \
    -siteMaxAlleleFreq 0.99 \
    -removeSitesWithIndels true \
    -endPlugin \
    -filterAlign \
    -subsetSites 300000 \
    -step \
    -export $VCF_OUTMERGED/subset4Kinship_hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr${CHROM}_imputed_sorted.vcf \
    -exportType VCF
  echo "I just finished chr ${CHROM}"
done


# Merge vcf files by chromosome
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
  -debug $VCF_OUTMERGED/debug_merge_vcf_nam_subset.log \
  -Xmx200g \
  -maxThreads 30 \
  -fork1 -importGuess $VCF_OUTMERGED/subset4Kinship_hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr1_imputed_sorted.vcf -noDepth \
  -fork2 -importGuess $VCF_OUTMERGED/subset4Kinship_hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr2_imputed_sorted.vcf -noDepth \
  -fork3 -importGuess $VCF_OUTMERGED/subset4Kinship_hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr3_imputed_sorted.vcf -noDepth \
  -fork4 -importGuess $VCF_OUTMERGED/subset4Kinship_hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr4_imputed_sorted.vcf -noDepth \
  -fork5 -importGuess $VCF_OUTMERGED/subset4Kinship_hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr5_imputed_sorted.vcf -noDepth \
  -fork6 -importGuess $VCF_OUTMERGED/subset4Kinship_hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr6_imputed_sorted.vcf -noDepth \
  -fork7 -importGuess $VCF_OUTMERGED/subset4Kinship_hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr7_imputed_sorted.vcf -noDepth \
  -fork8 -importGuess $VCF_OUTMERGED/subset4Kinship_hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr8_imputed_sorted.vcf -noDepth \
  -fork9 -importGuess $VCF_OUTMERGED/subset4Kinship_hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr9_imputed_sorted.vcf -noDepth \
  -fork10 -importGuess $VCF_OUTMERGED/subset4Kinship_hmp321_282_agpv5_inbred_hybrid_b73_phz51_phb47_uplifted_merged_chr10_imputed_sorted.vcf -noDepth \
  -combine11 -input1 -input2 -input3 -input4 -input5 -input6 -input7 -input8 -input9 -input10 \
  -mergeGenotypeTables \
  -export $VCF_OUTMERGED/subset4Kinship_hmp321_282_agpv5_inbred_hybrid_all_chroms.vcf \
  -exportType VCF



# all files backuped within: /data1/users/mbb262/te_ase_nam/genotypes