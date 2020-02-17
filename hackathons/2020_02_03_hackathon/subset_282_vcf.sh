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
for CHROM in {1..10}
do
  echo "Start data subsetting on chromosome ${CHROM}"
  bcftools index -c /workdir/mbb262/AGPv4/hmp321_282_agpv4_merged_chr${CHROM}.imputed.vcf.gz
  bcftools view -S /workdir/mbb262/goodman282/282set_lines.txt hmp321_282_agpv4_merged_chr${CHROM}.imputed.vcf.gz -Oz -o /workdir/mbb262/goodman282/hmp321_282_agpv4_merged_chr${CHROM}_imputed_goodman282.vcf.gz
  echo "End data subsetting on chromosome ${CHROM}"
done

# Filter SNPs based on imputation accuracy
for CHROM in {1..10}
do
  echo "Start data filtering on chromosome ${CHROM}"
  bcftools filter -e 'DR2<=0.8' hmp321_282_agpv4_merged_chr${CHROM}_imputed_goodman282.vcf.gz -o hmp321_282_agpv4_merged_chr${CHROM}_imputed_filteredDR2_goodman282.vcf.gz
  echo "End data filtering on chromosome ${CHROM}"
done

# Run SNP GWAS within tassel 4
for CHROM in {1..10}
do
  /programs/tassel4-standalone/run_pipeline.pl \
    -Xmx220g \
    -fork1 -importGuess /workdir/mbb262/goodman282/hmp321_282_agpv4_merged_chr${CHROM}_imputed_filteredDR2_goodman282.vcf.gz \
    -fork2 -r /workdir/mbb262/goodman282/vcf_intersect_mdp_traits.txt \
    -fork3 -ck -fork1 \
    -combine4 -input1 -input2 -intersect \
    -combine5 -input5 -input3 \
    -mlm -mlmVarCompEst P3D -mlmCompressionLevel None \
    -export /workdir/mbb262/goodman282/gwas_results/MLM_goodman282_kinship \
    -runfork1 -runfork2 -runfork3
    echo "End mlm on chromosome ${CHROM}"
done


# ------------------------------
# Run TASSEL 4 MLM on chr10 only
# ------------------------------

# Generate kinship matrix for NAM
/programs/tassel-5-standalone/run_pipeline.pl \
  -debug /workdir/mbb262/goodman282/debug_kinship_chr10.log \
  -Xmx220g \
  -importGuess /workdir/mbb262/goodman282/hmp321_282_agpv4_merged_chr10_imputed_filteredDR2_goodman282.vcf \
  -KinshipPlugin \
  -method Centered_IBS \
  -endPlugin \
  -export /workdir/mbb262/goodman282/goodman282_kinship_chr10.txt \
  -exportType SqrMatrix

# Run mlm in tassel4
/home/mbb262/bioinformatics/tassel4-standalone/run_pipeline.pl \
  -Xmx500g \
  -fork1 -importGuess /workdir/mbb262/goodman282/hmp321_282_agpv4_merged_chr10_imputed_filteredDR2_goodman282.vcf \
  -fork2 -importGuess /workdir/mbb262/goodman282/vcf_intersect_mdp_traits.txt \
  -fork3 -importGuess /workdir/mbb262/goodman282/goodman282_kinship_chr10.txt \
  -combine4 -input1 -input2 -intersect \
  -combine5 -input4 -input3 \
  -mlm -mlmVarCompEst P3D -mlmCompressionLevel None \
  -mlmOutputFile /workdir/mbb262/goodman282/gwas_results/MLM_goodman282_kinship \
  -runfork1 -runfork2 -runfork3



/workdir/mbb262/nam_founders/filtered_nam_founders








