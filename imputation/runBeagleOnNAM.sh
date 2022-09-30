#! /bin/bash

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


# Beagle imputation on vcf file
java \
  -Xmx450g \
  -jar /programs/beagle41/beagle41.jar \
  gt=/home/mbb262/tassel_gwa/hackathon_2019_10-24/ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023_main200_MAF035_MAF000001_NAMAmesintersect_noMissingAbove.36.vcf \
  ne=1000 \
  impute=true \
  nthreads=111 \
  out=/home/mbb262/tassel_gwa/hackathon_2019_10-24/intersectJoin_NAMbyAmes_shared_beagleImputed.vcf


# Beagle imputation on NAM founders vcf file
java \
  -Xmx450g \
  -jar /programs/beagle41/beagle41.jar \
  gt=/home/mbb262/tassel_gwa/data/genos/ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023_Nam_founders_sharedNAMAmesSites.vcf \
  ne=1000 \
  impute=true \
  nthreads=111 \
  out=/home/mbb262/tassel_gwa/data/genos/intersectJoin_NAMbyAmes_shared_beagleImputed_NAMfounders.vcf


# Beagle imputation on pre-imputed NAM families and non-imputed IBM
java \
  -Xmx450g \
  -jar /programs/beagle41/beagle41.jar \
  nthreads=60 \
  niterations=2 \
  gt=/home/mbb262/tassel_gwa/data/genos/intersectJoin_NAMbyAmes_shared40k_allNAM_IBM.vcf \
  ne=1000 \
  impute=true \
  out=/home/mbb262/tassel_gwa/data/genos/intersectJoin_NAMbyAmes_shared40k_allNAM_IBM_beagleImputation
