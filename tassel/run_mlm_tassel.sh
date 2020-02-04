#!/bin/bash

# ----------------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2019-11-26 
#
# Description 
# 	- General MLM script for TASSEL where y is 5-6 phenotypes with increasing 
#	  levels of complexity in NAM
#   - y ~ reduced nam SNP set + kinship matrix 
# 	- Model takes forever to run, kinship matrix function is fast
# ----------------------------------------------------------------------


# Generate kinship matrix for NAM
/programs/tassel-5-standalone/run_pipeline.pl \
  -debug /home/mbb262/tassel_gwa/results/2019_12_04_ames2nam_comparisons/debug_kinshipMatrix.log \
  -Xmx120g \
  -importGuess /home/mbb262/tassel_gwa/data/ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023_NAM_only_main200_maf035.vcf \
  -KinshipPlugin \
  -method Centered_IBS \
  -endPlugin \
  -export /home/mbb262/tassel_gwa/results/2019_12_04_ames2nam_comparisons/NAM_only_main200_maf035_kinshipMatrix.txt \
  -exportType SqrMatrix


# # Run MLM with kinship matrix for
# #	y ~ reduced nam SNP set + kinship matrix
/programs/tassel-5-standalone/run_pipeline.pl \
  -debug /home/mbb262/tassel_gwa/results/2019_12_04_ames2nam_comparisons/debug_MLM_NAMkinship_2019_12_04.log \
  -Xmx420g \
  -maxThreads 83 \
  -fork1 -importGuess /home/mbb262/tassel_gwa/data/ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023_NAM_only_main200_maf035.vcf \
  -fork2 -r /home/mbb262/tassel_gwa/data/range_complexity_phenos_wallace.txt \
  -fork3 -r /home/mbb262/tassel_gwa/data/ames2nam_3gPCs.txt \
  -fork4 -k /home/mbb262/tassel_gwa/results/2019_12_04_ames2nam_comparisons/NAM_only_main200_maf035_kinshipMatrix.txt \
  -combine5 -input1 -input2 -input3 -intersect \
  -combine6 -input5 -input4 \
  -mlm -mlmVarCompEst P3D -mlmCompressionLevel None \
  -export /home/mbb262/tassel_gwa/results/2019_12_04_ames2nam_comparisons/MLM_NAMkinship_2019_12_04 \
  -runfork1 -runfork2 -runfork3
