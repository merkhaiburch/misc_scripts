#!/bin/bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-06-15 
#
# Description 
#   - Run SKmer on subsampled PanAnd genomes and create a 
#   - phylogenetic tree
# ---------------------------------------------------------------

# Add paths to the key programs
export PATH=/programs/jellyfish-2.3.0/bin:$PATH
export PATH=/programs/mash-Linux64-v2.1:$PATH
export PATH=/programs/seqtk:$PATH

# test to see if the paths were added
jellyfish --version
mash --version
seqtk


# -------------------------------
# Run skmer on 100k random reads
# -------------------------------

# Create a reference 
skmer \
    reference \
    /workdir/hackathon/subsampled_seq2 \
    -p 20 \
    -l /workdir/hackathon/skmer_results/panand16results \
    -o /workdir/hackathon/skmer_results/panand_skmer_16_genomes

# Compute distances
skmer \
    distance \
    /workdir/hackathon/skmer_results/panand16results \
    -t \
    -o /workdir/hackathon/skmer_results/jc-dist-mat-pandand16


# -------------------------------
# Run skmer on 1M random reads
# -------------------------------

# Create a reference 
skmer \
    reference \
    /workdir/hackathon/subsampled_seq2_1m \
    -p 20 \
    -l /workdir/hackathon/skmer_results/panand16results_1m \
    -o /workdir/hackathon/skmer_results/panand_skmer_16_genomes_1m

# Compute distances
skmer \
    distance \
    /workdir/hackathon/skmer_results/panand16results_1m \
    -t \
    -o /workdir/hackathon/skmer_results/jc-dist-mat-pandand16_1m


# -------------------------------------------
# Run skmer on 1M random reads + new genomes
# -------------------------------------------

# Create a reference 
skmer \
    reference \
    /workdir/hackathon/seq3/one_mil \
    -p 20 \
    -l /workdir/hackathon/skmer_results/panand16results_1m_plus_new \
    -o /workdir/hackathon/skmer_results/panand_skmer_16_genomes_1m_plus_new

# Compute distances
skmer \
    distance \
    /workdir/hackathon/skmer_results/panand16results_1m_plus_new \
    -t \
    -p 20 \
    -o /workdir/hackathon/skmer_results/jc-dist-mat-pandand16_1m_plus_new

