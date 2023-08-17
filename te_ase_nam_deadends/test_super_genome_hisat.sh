#!/usr/bin/env bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-04-27 
# Updated... 2022-07-22
#
# Description:
# - Test script to run minimap2 to align reads back to their super 
# 	founder genome using minimap2
# ---------------------------------------------------------------


# ------------------------------------------------------------------------------------------
# Align RNAseq reads with minimap2 to transcriptome
# ------------------------------------------------------------------------------------------

# Gather genomes -------------------------------------------------------------------------------

# make directories, gather transcripts
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/references
cd /workdir/mbb262/nam_hybrid_rnaseq/references
scp -r mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/te_ase_nam/b73_trans_aligned_to_nam_aimee/renamed_alignments ./

# Gather test rnaseq reads -------------------------------------------------------------------------

# All growing point tissues
# B73/ky21: MS21R233
# B73: MS21R001
# Ky21: MS21R021

# Make directories
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/reads/raw
cd /workdir/mbb262/nam_hybrid_rnaseq/reads/raw

# B73/ky21: MS21R233
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R233/MS21R233_CKDL210018333-2a-AK17213-7UDI233_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R233/MS21R233_CKDL210018333-2a-AK17213-7UDI233_HH5V7DSX2_L1_2.fq.gz

# B73: MS21R020
# Note B73: MS21R001 --> failed library I think
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R020/MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R020/MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_2.fq.gz

# Ky21: MS21R021
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R021/MS21R021_CKDL210018333-2a-AK2147-7UDI223_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R021/MS21R021_CKDL210018333-2a-AK2147-7UDI223_HH5V7DSX2_L1_2.fq.gz

# unit test reads
cp /home/mbb262/git_projects/te_ase_nam/data/unit_tests* ./
gzip unit_tests_*


# Run TrimGalore ---------------------------------------------------------------------------------

# Make directories
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/reads/raw
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/reads/trimmed
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/reads/trimming_reports
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/output/fastqc

## Parameters ----
N_THREADS=20
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
FASTQC_RAW_DIR=$PROJ_DIR/reads/raw 
FASTQC_OUT=$PROJ_DIR/reads/trimming_reports
FASTQC_TRIM=$PROJ_DIR/reads/trimmed/


## Trim Galore! for paired end reads -------------------------------------------------------------

# For all files in directory, including unit tests
for i in $FASTQC_RAW_DIR/*_1.fq.gz
do
    SAMPLE=$(basename ${i} _1.fq.gz)

    echo "Trimming: " ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz

    trim_galore \
    --cores $N_THREADS \
    --paired $FASTQC_RAW_DIR/${SAMPLE}_1.fq.gz $FASTQC_RAW_DIR/${SAMPLE}_2.fq.gz \
    --length 148 \
    --output_dir $FASTQC_TRIM \
    --quality 20 \
    --fastqc \
    --fastqc_args "-o $FASTQC_OUT" \
    --basename $SAMPLE
done

# move summary outputs to a differnt directory
mkdir /workdir/mbb262/nam_hybrid_rnaseq/reads/trimming_reports
mv /workdir/mbb262/nam_hybrid_rnaseq/reads/trimmed/*_trimming_report.txt /workdir/mbb262/nam_hybrid_rnaseq/reads/trimming_reports


## Parse summary output to table ---------------------------------------------------------------
export LC_ALL=en_US.UTF-8
export PYTHONPATH=/programs/multiqc-1.10.1/lib64/python3.6/site-packages:/programs/multiqc-1.10.1/lib/python3.6/site-packages
export PATH=/programs/multiqc-1.10.1/bin:$PATH

multiqc /workdir/mbb262/nam_hybrid_rnaseq/reads/trimming_reports
mv /workdir/mbb262/nam_hybrid_rnaseq/reads/trimmed/multiqc_* /workdir/mbb262/nam_hybrid_rnaseq/output/

# Go from sam to bam ---------------------------------------------------------------------------
cd /workdir/mbb262/nam_hybrid_rnaseq/references/renamed_alignments
samtools view -bS Zm-B73-REFERENCE-NAM-5.0_B73v5_mRNA_axSplice_eqx.sam > Zm-B73-REFERENCE-NAM-5.0_B73v5_mRNA_axSplice_eqx.sam.bam
samtools view -bS Zm-Ky21-REFERENCE-NAM-1.0_B73v5_mRNA_axSplice_eqx.sam > Zm-Ky21-REFERENCE-NAM-1.0_B73v5_mRNA_axSplice_eqx.sam.bam


# bam to fa -------------------------------------------------------------------------------------
samtools fasta Zm-B73-REFERENCE-NAM-5.0_B73v5_mRNA_axSplice_eqx.sam.bam > b73_b73.fa
samtools fasta Zm-Ky21-REFERENCE-NAM-1.0_B73v5_mRNA_axSplice_eqx.sam.bam > b73_ky21.fa


# Gather only canonical transcripts from B73 ----------------------------------------------------
# to be replaced with better gene choices later
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts

# subsample fasta files to just canonical transcripts 
samtools faidx b73_b73.fa -r Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts -o b73_b73_canonical.fa 
samtools faidx b73_ky21.fa -r Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts -o ky21_b73_canonical.fa 

# Combine b73 and ky21 fastas once appending on helpful names
sed '/^>/s/$/_b73/' b73_b73_canonical.fa > b73_b73_canonical_named.fa 
sed '/^>/s/$/_ky21/' ky21_b73_canonical.fa > ky21_b73_canonical_named.fa
cat b73_b73_canonical_named.fa ky21_b73_canonical_named.fa > b73_ky21_combined_canonical_trans.fa


## Run minimap2 -sr for rnaseq reads onto transcriptomes --------------------------------------------

mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/output/minimap_alignments

# variables
N_THREADS=20
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
ALIGN_OUT=$PROJ_DIR/output/minimap_alignments
REF_TRANS=$PROJ_DIR/references/renamed_alignments
FASTQC_TRIM=$PROJ_DIR/reads/trimmed

# unit tests to b73 #1
/programs/minimap2-2.17/minimap2 -ax sr \
    -t $N_THREADS \
    --max-qlen 350 \
    $REF_TRANS/b73_b73_canonical_named.fa \
    $FASTQC_TRIM/unit_tests_modified_val_1.fq.gz \
    $FASTQC_TRIM/unit_tests_modified_val_2.fq.gz \
    > $ALIGN_OUT/minimap_sr_unitTest_to_b73_canonical_all.sam

# unit tests to b73 #2
/programs/minimap2-2.17/minimap2 -ax sr \
    -t $N_THREADS \
    --max-qlen 350 \
    $REF_TRANS/b73_b73_canonical_named.fa \
    $FASTQC_TRIM/unit_tests_count_val_1.fq.gz \
    $FASTQC_TRIM/unit_tests_count_val_2.fq.gz \
    > $ALIGN_OUT/minimap_sr_unitTest_to_b73_canonical_count.sam

# unit tests to b73 #3
/programs/minimap2-2.17/minimap2 -ax sr \
    -t $N_THREADS \
    --max-qlen 350 \
    $REF_TRANS/b73_b73_canonical_named.fa \
    $FASTQC_TRIM/unit_tests_nocount_val_1.fq.gz \
    $FASTQC_TRIM/unit_tests_nocount_val_2.fq.gz \
    > $ALIGN_OUT/minimap_sr_unitTest_to_b73_canonical_nocount.sam


# Make a loop that goes through all reads 





# Filter or uniquely count though sam files that:
# keep things with less than 3 mismatches



