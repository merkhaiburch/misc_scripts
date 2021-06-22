#!/bin/bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-06-09 
#
# Description 
#   - Calculate a genomic relationship matrix for NAM and 282
#   - Use PHG imputed SNPs for:
#   - 1) genome wide SNPs 2) SNPs only within specific intervals
# ---------------------------------------------------------------

# --------------------
#  NAM SNPs and GRM
# --------------------

# Filter data using tassel using MAF then randomly subset 10M sites genome wide
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/debug_calculate_grm_nam.log \
    -Xmx1000g \
    -maxThreads 86 \
    -importGuess /workdir/mbb262/NAM_phg_snps.vcf -noDepth \
    -FilterSiteBuilderPlugin \
    -siteMinAlleleFreq 0.01 \
    -siteMaxAlleleFreq 0.99 \
    -removeSitesWithIndels true \
    -endPlugin \
    -filterAlign \
    -subsetSites 10000000 \
    -step \
    -export /workdir/mbb262/filtered_NAM_phg_snps.vcf \
    -exportType VCF

gzip /workdir/mbb262/filtered_NAM_phg_snps.vcf

# calculate kinship genome wide using a subsetted vcf
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/debug_calculate_grm_kinship_nam.log \
    -Xmx1000g \
    -maxThreads 86 \
    -importGuess /workdir/mbb262/filtered_NAM_phg_snps.vcf.gz -noDepth \
    -KinshipPlugin \
    -method Centered_IBS \
    -endPlugin \
    -export /workdir/mbb262/kinship_filtered_NAM_phg_snps.txt

# if the above doesn't work, try a kinship matirx on the whole genome file


# -----------------------------------------
#  NAM SNPs and GRM --> HARE regions only
# -----------------------------------------

# # bed file from HARE genes in v3, filtering SNPs within these intervals
# # Filter data using tassel using MAF 
# /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
#     -debug /workdir/mbb262/nam_hare_downsample_debug.log \
#     -Xmx1000g \
#     -maxThreads 86 \
#     -importGuess /workdir/workdir/NAM_phg_snps.vcf -noDepth \
#     -FilterSiteBuilderPlugin \
#     -siteMinAlleleFreq 0.01 \
#     -siteMaxAlleleFreq 0.99 \
#     -removeSitesWithIndels true \
#     -bedFile ./hare_v5_intervals.bed \
#     -endPlugin \
#     -export /workdir/mbb262/filtered_hare_nam_phg_snps.vcf \
#     -exportType VCF

# gzip /workdir/mbb262/filtered_hare_nam_phg_snps.vcf

# calculate kinship genome wide using a subsetted vcf
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/nam_hare_calc_grm_debug.log \
    -Xmx1000g \
    -maxThreads 86 \
    -importGuess /workdir/mbb262/filtered_NAM_phg_snps.vcf.gz -noDepth \
    -KinshipPlugin \
    -method Centered_IBS \
    -endPlugin \
    -export /workdir/mbb262/kinship_filtered_nam_hare_phg_snps.txt


# --------------------
#  282 SNPs and GRM
# --------------------

# Filter data using tassel using MAF 
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/goodman_downsample_debug.log \
    -Xmx500g \
    -maxThreads 60 \
    -importGuess /workdir/ag2484/goodman_phg_snps.txt.vcf -noDepth \
    -filterAlign \
    -filterAlignMinFreq 0.01 \
    -filterAlignMaxFreq 0.99 \
    -export /workdir/mbb262/filtered_goodman_phg_snps.vcf.gz \
    -exportType VCF


# calculate kinship genome wide using a subsetted vcf
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/goodman_calc_grm_debug.log \
    -Xmx500g \
    -maxThreads 60 \
    -importGuess /workdir/mbb262/filtered_goodman_phg_snps.vcf.gz -noDepth \
    -KinshipPlugin \
    -method Centered_IBS \
    -endPlugin \
    -export /workdir/mbb262/kinship_filtered_goodman_phg_snps.txt


# -----------------------------------------
#  282 SNPs and GRM --> HARE regions only
# -----------------------------------------

# Filter data using tassel using MAF 
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/goodman_hare_downsample_debug.log \
    -Xmx500g \
    -maxThreads 60 \
    -importGuess /workdir/ag2484/goodman_phg_snps.txt.vcf -noDepth \
    -FilterSiteBuilderPlugin \
    -siteMinAlleleFreq 0.01 \
    -siteMaxAlleleFreq 0.99 \
    -removeSitesWithIndels true \
    -bedFile hare_v5_intervals.bed \
    -endPlugin \
    -export /workdir/mbb262/filtered_hare_goodman_phg_snps.vcf \
    -exportType VCF

gzip /workdir/mbb262/filtered_hare_goodman_phg_snps.vcf

# calculate kinship genome wide using a subsetted vcf
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/goodman_hare_calc_grm_debug.log \
    -Xmx500g \
    -maxThreads 60 \
    -importGuess /workdir/mbb262/filtered_hare_goodman_phg_snps.vcf.gz -noDepth \
    -KinshipPlugin \
    -method Centered_IBS \
    -endPlugin \
    -export /workdir/mbb262/kinship_filtered_goodman_hare_phg_snps.txt








# Split results into different chromosomes
bgzip -c /workdir/ag2484/NAM_phg_snps.vcf > /workdir/mbb262/NAM_phg_snps.vcf.gz
tabix -p vcf /workdir/mbb262/NAM_phg_snps.vcf.gz

# loop extract's each chromosome's info
for CHROM in {1..10}
do
    echo "Start analyzing chromosome ${CHROM}"
    tabix /workdir/mbb262/NAM_phg_snps.vcf.gz chr${CHROM} > nam_chr${CHROM}.vcf
    echo "End analyzing chromosome ${CHROM}"
done




# Other non-functional code
# randomly subset 10M sites genome wide using jvarkit/downsamplevcf
java -jar /home/mbb262/bioinformatics/jvarkit/dist/downsamplevcf.jar \
    /workdir/mbb262/filtered_goodman_phg_snps.vcf \
    -o downsampled_10M_filtered_goodman_phg_snps.vcf.gz \
    -n 10000000

bcftools view --header-only /workdir/mbb262/filtered_goodman_phg_snps.vcf > subsample_goodman.vcf
bcftools view --no-header input.vcf | awk '{printf("%f\t%s\n",rand(),$0);}' | sort -t $'\t' -k1,1g | cut -f2-  | head -n 1000 >>  subsample.vcf
