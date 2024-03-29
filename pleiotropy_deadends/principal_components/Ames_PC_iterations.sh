#!/bin/bash

# -----------------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-09-04
#
# Script to run TASSEL through the command-line
# to calculate a distance/relationship matrix
# for the Ames population
# For local windows try # 10Mb per window
# -----------------------------------------------


# -------------------------------------
#        Get clean data
# -------------------------------------


# Import unimputed Ames data
# Filter on min count
/programs/tassel-5-standalone/run_pipeline.pl \
  -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_filteredGenotypes_log.txt \
  -Xmx250g \
  -importGuess /home/panzea/processedData/genotypes/ZeaGBSv27/ZeaGBSv27_20171204_AGPv4_Ames.vcf \
  -FilterSiteBuilderPlugin \
  -siteMinCount 2500 \
  -siteMinAlleleFreq 0.0001 \
  -endPlugin \
  -FilterTaxaBuilderPlugin \
  -minNotMissing 0.55 \
  -endPlugin \
  -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes.vcf -exportType VCF

# Commands to check depth in vcf files
grep "DP=0" ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes.vcf > lala.txt
wc -l lala.txt
less -S ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes1.vcf


# ------------------------------------------
#       No imputation pipeline
# ------------------------------------------

# PCs directly on unimputed VCF
/programs/tassel-5-standalone/run_pipeline.pl \
 -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_PCsonVCF_no_imputation_log.txt \
 -Xmx120g \
 -importGuess /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes.vcf \
 -PrincipalComponentsPlugin \
 -covariance true \
 -endPlugin \
 -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/no_imputation_output

# Calculate distance matrix on unimputed data
/programs/tassel-5-standalone/run_pipeline.pl \
  -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_matricesForPCsDistance_no_imputation_log.txt \
  -Xmx120g \
  -importGuess /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes.vcf \
  -DistanceMatrixPlugin \
  -endPlugin \
  -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_no_imputation_genotypeDistance.txt

# Calculate kinship matrix on unimputed data
/programs/tassel-5-standalone/run_pipeline.pl \
  -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_matricesForPCsKinship_no_imputation_log.txt \
  -Xmx120g \
  -importGuess /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes.vcf \
  -KinshipPlugin \
  -method Centered_IBS \
  -endPlugin \
  -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_no_imputation_kinshipMatrix.txt \
  -exportType SqrMatrix

# Remove headers generated by tassel
# Example:
##IBS_Distance_Matrix.AverageTotalSites=166219.08371645436
##IBS_Distance_Matrix.NumAlleles=3
##IBS_Distance_Matrix.TrueIBS=false
##Matrix_Type=IBS_Distance_Matrix
# 3545
# grep -vwE "##IBS_Distance_Matrix" ZeaGBSv27_20171204_AGPv4_Ames_genotype_distance.txt > ZeaGBSv27_20171204_AGPv4_Ames_genotype_distance2.txt
# grep -vwE "##Matrix_Type" ZeaGBSv27_20171204_AGPv4_Ames_genotype_distance2.txt > ZeaGBSv27_20171204_AGPv4_Ames_genotype_distance3.txt
# grep -vwE "3545" ZeaGBSv27_20171204_AGPv4_Ames_genotype_distance3.txt > ZeaGBSv27_20171204_AGPv4_Ames_genotype_distance4.txt


# ----------------------------------
#     Mean Imputation Pipeline
# ----------------------------------

# Do imputation (mean)
# Export mean imputed vcf file
/programs/tassel-5-standalone/run_pipeline.pl \
  -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_filteredGenotypesImputed_log.txt \
  -Xmx120g \
  -importGuess /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes.vcf \
  -ImputationPlugin \
  -ByMean true \
  -endPlugin \
  -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypesImputed -exportType ReferenceProbability


# PCs directly on VCF
/programs/tassel-5-standalone/run_pipeline.pl \
 -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_PCsonVCF_log.txt \
 -Xmx120g \
 -importGuess /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypesImputed.txt \
 -PrincipalComponentsPlugin \
 -covariance true \
 -endPlugin \
 -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/impute_output

# Cannot calculate a distance or kinship matrix on a numerical imputation output file


# -----------------------------
#  Beagle imputation pipeline
# -----------------------------

# Beagle imputation on vcf file
java \
  -Xmx500g \
  -jar /programs/beagle41/beagle41.jar \
  gt=/home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes.vcf \
  ne=1000 \
  impute=true \
  nthreads=95 \
  out=/home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes_beagleImputation

# PCs directly on VCF
/programs/tassel-5-standalone/run_pipeline.pl \
 -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_PCsonVCF_beagle_log.txt \
 -Xmx120g \
 -importGuess /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/beagle_imputed_ames_snps.vcf \
 -PrincipalComponentsPlugin \
 -covariance true \
 -endPlugin \
 -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/beagle_impute_output \

# Calculate distance matrix on beagle imputed data
/programs/tassel-5-standalone/run_pipeline.pl \
  -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_matricesForPCsDistance_beagle_log.txt \
  -Xmx120g \
  -importGuess /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/beagle_imputed_ames_snps.vcf \
  -DistanceMatrixPlugin \
  -endPlugin \
  -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_beagle_imputation_genotypeDistance.txt

# Calculate kinship matrix on beagle imputed data
/programs/tassel-5-standalone/run_pipeline.pl \
  -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_matricesForPCsKinship_beagle_log.txt \
  -Xmx120g \
  -importGuess /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/beagle_imputed_ames_snps.vcf \
  -KinshipPlugin \
  -method Centered_IBS \
  -endPlugin \
  -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_beagle_imputation_kinshipMatrix.txt \
  -exportType SqrMatrix

# calculate PCs on distance matrix and kinship matrix in R


# ---------------------------
#   KNNI imputation pipeline
# ---------------------------

# Did imputation in tassel GUI
# Exported as:

# PCs directly on VCF
/programs/tassel-5-standalone/run_pipeline.pl \
 -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_PCsonVCF_knni_log.txt \
 -Xmx120g \
 -importGuess /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/fillin_imputed_ames_sep20.hmp.txt \
 -PrincipalComponentsPlugin \
 -covariance true \
 -endPlugin \
 -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/knni_impute_output \

# Calculate distance matrix on KNNI imputed data
/programs/tassel-5-standalone/run_pipeline.pl \
  -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_matricesForPCsDistance_knni_log.txt \
  -Xmx120g \
  -importGuess /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/fillin_imputed_ames_sep20.hmp.txt \
  -DistanceMatrixPlugin \
  -endPlugin \
  -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_knni_imputation_genotypeDistance.txt

# Calculate kinship matrix on KNNI imputed data
/programs/tassel-5-standalone/run_pipeline.pl \
  -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_matricesForPCsKinship_knni_log.txt \
  -Xmx120g \
  -importGuess /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/fillin_imputed_ames_sep20.hmp.txt \
  -KinshipPlugin \
  -method Centered_IBS \
  -endPlugin \
  -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_knni_imputation_kinshipMatrix.txt \
  -exportType SqrMatrix


# ---------------------------
#   FILLIN imputation pipeline
# ---------------------------

# Made haplotypes in tassel GUI
# Did imputation in tassel GUI
# Exported file as: 

# Make haplotypes with KNNI


# PCs directly on VCF
/programs/tassel-5-standalone/run_pipeline.pl \
 -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_PCsonVCF_fillin_log.txt \
 -Xmx120g \
 -importGuess /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/fillin_imputed_ames_sep20.hmp.txt \
 -PrincipalComponentsPlugin \
 -covariance true \
 -endPlugin \
 -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/fillin_impute_output \

# Calculate distance matrix on fillin imputed data
/programs/tassel-5-standalone/run_pipeline.pl \
  -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_matricesForPCsDistance_fillin_log.txt \
  -Xmx120g \
  -importGuess /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/fillin_imputed_ames_sep20.hmp.txt \
  -DistanceMatrixPlugin \
  -endPlugin \
  -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_fillin_imputation_genotypeDistance.txt

# Calculate kinship matrix on fillin imputed data
/programs/tassel-5-standalone/run_pipeline.pl \
  -debug /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/Ames_matricesForPCsKinship_fillin_log.txt \
  -Xmx120g \
  -importGuess /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/fillin_imputed_ames_sep20.hmp.txt \
  -KinshipPlugin \
  -method Centered_IBS \
  -endPlugin \
  -export /home/mbb262/tassel_gwa/ames_FA_tests/Sep19_whole_pipeline/ZeaGBSv27_20171204_AGPv4_Ames_fillin_imputation_kinshipMatrix.txt \
  -exportType SqrMatrix


# Manually deleted headers in all files before importing into R


