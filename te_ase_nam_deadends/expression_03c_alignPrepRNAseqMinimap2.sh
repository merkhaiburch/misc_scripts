#!/usr/bin/env bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-09-12 
# Updated... 2023-02-22
#
# Description: 
# Last script ran minimap2 to align B73 canonical transcripts back 
# to each individual NAM inbred founder
# AKA: "fish" out the same B73 transcript in each reference genome 
#
# This script formats those alignments by:
# - picking the primary alignment, and finding its underlying sequence
# 
# Next script is created by kotlin and pairs each sample with its
# corresponding reference transcriptome and generates a minimap script to
# align such reads to these formatted reference transcriptomes
# ---------------------------------------------------------------

# Set project directory
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq


# Make directories ----------------------------------------------
mkdir -p $PROJ_DIR/references
mkdir -p $PROJ_DIR/references/bam_transcriptomes
mkdir -p $PROJ_DIR/references/fa_transcriptomes
mkdir -p $PROJ_DIR/references/nam_aligned_transcriptomes
mkdir -p $PROJ_DIR/references/bam_transcriptomes/primary


# Make directories, gather transcriptome alignments -------------

cd $PROJ_DIR/references
scp -r mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/te_ase_nam/nam_transcriptomes_primary/* $PROJ_DIR/references/nam_aligned_transcriptomes
ll $PROJ_DIR/references/nam_aligned_transcriptomes


# SAM to BAM to FA to subsampled FA transcriptomes --------------

# Variables for converting between file types
SAM_TRANSCRIPTOME_DIR=$PROJ_DIR/references/nam_aligned_transcriptomes
BAM_TRANSCRIPTOME_DIR=$PROJ_DIR/references/bam_transcriptomes
FA_TRANSCRIPTOME_DIR=$PROJ_DIR/references/fa_transcriptomes
PRIMARY_ALIGN_DIR=$PROJ_DIR/references/bam_transcriptomes/primary

# Loop through all transcriptomes and format them
for i in $SAM_TRANSCRIPTOME_DIR/*.sam
do
    SAMPLE=$(basename ${i} .sam)

    # echo "SAM to BAM to FA to subsampled FA: " $SAM_TRANSCRIPTOME_DIR/${SAMPLE}.sam
    echo -e "SAM to BAM to FA to FA subset: ${SAMPLE}"

    # sam to bam
    /programs/samtools-1.15.1-r/bin/samtools view -bS $SAM_TRANSCRIPTOME_DIR/${SAMPLE}.sam > $BAM_TRANSCRIPTOME_DIR/${SAMPLE}.bam

    # Index bam
    /programs/samtools-1.15.1-r/bin/samtools index $BAM_TRANSCRIPTOME_DIR/${SAMPLE}.bam

    # Get primary alignments from bam
    /programs/samtools-1.15.1-r/bin/samtools view -e 'flag <= 16' $BAM_TRANSCRIPTOME_DIR/${SAMPLE}.bam -o $PRIMARY_ALIGN_DIR/primary_${SAMPLE}.bam

    # bam to bed
    bedtools bamtobed -i $PRIMARY_ALIGN_DIR/primary_${SAMPLE}.bam > $PRIMARY_ALIGN_DIR/primary_${SAMPLE}.bed
done


# bed to fasta
cd $PRIMARY_DIR
gunzip *0.fa.gz
find $PRIMARY_DIR -maxdepth 1 -name '*.fa' > $PRIMARY_ALIGN_DIR/refGenome.list
find $PRIMARY_ALIGN_DIR/ -maxdepth 1 -name '*.bed' > $PRIMARY_ALIGN_DIR/genomeBeds.list
cd $PRIMARY_ALIGN_DIR
parallel --link "bedtools getfasta -name -s -fi {1} -bed {2/} -fo {1/.}.fasta" :::: $PRIMARY_ALIGN_DIR/refGenome.list :::: $PRIMARY_ALIGN_DIR/genomeBeds.list 


# Do last bit of formatting to read names
for i in $PRIMARY_ALIGN_DIR/*.fasta
do
    SAMPLE=$(basename ${i} .fasta)

    # Remove where in each reference the seqeunce was pulled from
    sed 's/:.*//' $PRIMARY_ALIGN_DIR/${SAMPLE}.fasta > $PRIMARY_ALIGN_DIR/${SAMPLE}_shortNames.fasta
done






# # Clean up fa_transcriptomes directory
# mkdir /workdir/mbb262/nam_hybrid_rnaseq/references/fai_files
# mkdir /workdir/mbb262/nam_hybrid_rnaseq/references/unnamed_transcripts

# mv /workdir/mbb262/nam_hybrid_rnaseq/references/fa_transcriptomes/*fai /workdir/mbb262/nam_hybrid_rnaseq/references/fai_files
# mv /workdir/mbb262/nam_hybrid_rnaseq/references/fa_transcriptomes/*_eqx.fa /workdir/mbb262/nam_hybrid_rnaseq/references/unnamed_transcripts
# mv /workdir/mbb262/nam_hybrid_rnaseq/references/fa_transcriptomes/*_eqx_canonical.fa /workdir/mbb262/nam_hybrid_rnaseq/references/unnamed_transcripts




# # Align back to genomes ------------------------------------------

# # Create directories
# mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/output/minimap_alignments

# # Mapping parameters ---------------
# N_THREADS=45
# PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
# ALIGN_OUT=$PROJ_DIR/output/minimap_alignments
# REF_TRANS=$PROJ_DIR/references/combined_fa_transcriptomes
# FASTQ_TRIM_MERGE_DIR=$PROJ_DIR/reads/merged

# minimap2 -ax sr \
#     -t $N_THREADS \
#     $REF_TRANS/${inbred1} \
#     $FASTQ_TRIM_MERGE_DIR/${sample} > $ALIGN_OUT/${sample}_minimap2.sam")







