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

# Create a reference 
/home/mbb262/bioinformatics/Skmer/skmer \
    reference \
    /home/hackathon/subsampled_seq2 \
    -p 20 \
    panand_skmer_16_genomes \
    -l /workdir/hackathon/skmer_results \
    -o panand_skmer_16_genomes

# Compute distances
/home/mbb262/bioinformatics/Skmer/skmer \
    distance \
    panand_skmer_16_genomes \
    -t \
    -o jc-dist-mat-pandand16

# Run the query
/home/mbb262/bioinformatics/Skmer/skmer \
    query \
    qry.fastq \
    library