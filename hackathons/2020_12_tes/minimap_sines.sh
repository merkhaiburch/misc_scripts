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


# ----------------
# Minimap stuff
# ----------------

# Get help
# /programs/minimap2-2.17/minimap2 --help

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


# SINE finder parameters parameters
minimap2 -ax sr \
    -t 30 \
    -A 1 \
    -B 5 \
    -O 39,81 \
    -N 15000 \
    -f 0 \
    -k 4 \
    -w 1 \
    -n 2 \
    --secondary=yes \
    ./b73_v5.mmi \
    ./maize_TE_db_exemplars.sine.fa > sine_finder_align.sam


wc -l sine_finder_align.sam




