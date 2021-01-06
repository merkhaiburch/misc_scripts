#!/bin/bash

# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-12-15 
#
# Description 
#   - December 2020 hackathon: TE annotation group
#   - Minimap only known SINE TEs to maize v5, compare with sine-finder results
# ---------------------------------------------------------------

# Get B73 v5 
# wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
# gunzip *.gz


# ------------------
# Minimap2 on SINEs
# ------------------

# Set path to program
export PATH=/programs/minimap2-2.17:$PATH

# Index the genome to save time
minimap2 -d \
    ./b73_v5.mmi \
    ./Zm-B73-REFERENCE-NAM-5.0.fa 

# Using general usage parameters
minimap2 -ax sr \
    -t 30 \
    -A 1 \
    -B 5 \
    -O 39,81 \
    -N 15000 \
    -f 600000 \
    -k 4 \
    --secondary=yes \
    ./b73_v5.mmi \
    ./maize_TE_db_exemplars.sine.fa > sine_alignment.sam


# Count the number of matches
wc -l sine_alignment.sam

# convert sam to bam and load into r
samtools view -S -b \
    sine_alignment.sam > sine_alignment.bam
samtools view sine_alignment.bam | head

# Download v1 of maize genome
wget http://ftp.maizesequence.org/release-4a.53/assembly/ZmB73_AGPv1_genome.fasta.gz
gunzip *

# Index the genome
minimap2 -d \
    ./b73_v1.mmi \
    ./ZmB73_AGPv1_genome.fasta

# Run minimap on v1 of maize genome
minimap2 -ax sr \
    -t 30 \
    -A 1 \
    -B 5 \
    -O 39,81 \
    -N 15000 \
    -f 600000 \
    -k 4 \
    --secondary=yes \
    ./b73_v1.mmi \
    ./maize_TE_db_exemplars.sine.fa > v1_sine_alignment.sam

# convert sam to bam and load into r
samtools view -S -b \
    v1_sine_alignment.sam > v1_sine_alignment.bam
samtools view v1_sine_alignment.bam | head


# ------------------
# Minimap2 on Ds1
# ------------------


# Run minimap2 on Ds1 elements in v5
minimap2 -ax sr \
    -t 30 \
    -A 1 \
    -B 5 \
    -O 39,81 \
    -N 15000 \
    -f 600000 \
    -k 4 \
    --secondary=yes \
    ./b73_v5.mmi \
    ./ds1.fa > ds1_alignment.sam


# Count the number of matches
wc -l ds1_alignment.sam

# convert sam to bam and load into r
samtools view -S -b \
    ds1_alignment.sam > ds1_alignment.bam
samtools view ds1_alignment.bam | head


# Run minimap2 on Ds1 elements in v2
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-7/fasta/zea_mays/dna/Zea_mays.AGPv2.dna.toplevel.fa

# Index the genome
minimap2 -d \
    ./b73_v2.mmi \
    ./Zea_mays.AGPv2.dna.toplevel.fa

# -f 0 filter nothing
# -w 1 keep all matching kmers

minimap2 -ax sr \
    -t 30 \
    -A 1 \
    -B 5 \
    -O 39,81 \
    -N 15000 \
    -f 600000 \
    -k 4 \
    --secondary=yes \
    ./b73_v2.mmi \
    ./ds1.fa > v2_ds1_alignment.sam


# Count the number of matches
wc -l v2_ds1_alignment.sam

# convert sam to bam and load into r
samtools view -S -b \
    v2_ds1_alignment.sam > v2_ds1_alignment.bam
samtools view v2_ds1_alignment.bam | head

# ds1 to v1
minimap2 -ax sr \
    -t 30 \
    -A 1 \
    -B 5 \
    -O 39,81 \
    -N 15000 \
    -f 600000 \
    -k 4 \
    --secondary=yes \
    /workdir/mbb262/sine_results/b73_v1.mmi \
    ./ds1.fa > v1_ds1_alignment.sam


# Count the number of matches
wc -l v1_ds1_alignment.sam

# convert sam to bam and load into r
samtools view -S -b \
    v1_ds1_alignment.sam > v1_ds1_alignment.bam
samtools view v1_ds1_alignment.bam | head

