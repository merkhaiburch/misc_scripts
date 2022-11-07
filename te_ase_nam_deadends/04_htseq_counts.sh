#!/usr/bin/env bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-01-19 
#
# Description:
# Quantify counts of RNAseq reads using HTseq from 
#	hisat2 alignments
# ---------------------------------------------------------------

## Parameters ---------
N_THREADS=35
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
SAM_DIR=$PROJ_DIR/output/hisat2_alignments
GTF_FILE=$PROJ_DIR/references/b73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gtf


## Set environment ----
export PYTHONPATH=/programs/HTSeq-0.11.2/lib64/python3.6/site-packages/
export PATH=/programs/HTSeq-0.11.2/bin:$PATH


## Run HTSeq ----------
ls $PROJ_DIR/output/hisat2_alignments/*_val_output.sam | \
    parallel -j $N_THREADS  \
    "htseq-count \
    -f sam \
    -s yes \
    -m union \
    --nonunique all \
    {} \
    $GTF_FILE > {.}.counts"


## Move all counts to their own directory
mv $PROJ_DIR/output/hisat2_alignments/*.counts $PROJ_DIR/output/htseq_counts


# Gzip all sam files ---
# gzip $PROJ_DIR/output/hisat2_alignments/*.sam
