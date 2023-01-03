#!/usr/bin/env bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-11-12 
#
# Description 
#   - Convert SAM files to BAM files using samtools
# ---------------------------------------------------------------

## Parameters ----
N_THREADS=20
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq


## Export paths ----
export PATH=/programs/parallel/bin:$PATH


# === PARALLEL methods ==============================================

ls $PROJ_DIR/output/hisat2_alignments/*.sam.gz | \
    parallel -j $N_THREADS \
    "samtools view -S -b {} > {.}.bam"
