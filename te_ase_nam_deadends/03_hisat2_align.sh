#!/usr/bin/env bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-11-12 
# Updated... 2022-04-28
#
# Description 
#   - Run hisat2 to align reads back to founder genome and B73 v5
# ---------------------------------------------------------------

## Make directories ----------------
# mkdir /workdir/mbb262/nam_hybrid_rnaseq/hisat2_index
# mkdir /workdir/mbb262/nam_hybrid_rnaseq/output/hisat2_alignments
# mkdir /workdir/mbb262/nam_hybrid_rnaseq/output/summaries

## Gather paths --------------------
source /programs/HISAT2/hisat2.sh


## Index parameters ----------------
N_THREADS=38
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
GENOME_DIR=$PROJ_DIR/hisat2_index
GENOME_FA=$PROJ_DIR/references/b73/Zm-B73-REFERENCE-NAM-5.0.fa


## Build index ---------------------
hisat2-build \
    $GENOME_FA \
    $GENOME_DIR/b73


# Mapping parameters ---------------
N_THREADS=35
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
FASTQC_TRIM_DIR=$PROJ_DIR/reads/trimmed
GENOME_DIR=$PROJ_DIR/hisat2_index
GENOME_FA=$PROJ_DIR/references/Zm-B73-REFERENCE-NAM-5.0.fa
OUT_PATH=$PROJ_DIR/output/hisat2_alignments


## Run hisat2 ----------------------

# Gather base name for paired end read, 
# align both reads using hisat2 to b73
for i in $FASTQC_TRIM_DIR/*_1.fq.gz
do
    SAMPLE=$(basename ${i} _1.fq.gz)

    echo "Aligning: " ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz

    hisat2 \
        -p $N_THREADS \
        -x $GENOME_DIR/b73 \
        -1 $FASTQC_TRIM_DIR/${SAMPLE}_1.fq.gz \
        -2 $FASTQC_TRIM_DIR/${SAMPLE}_2.fq.gz \
        -S $OUT_PATH/${SAMPLE}_output.sam --new-summary \
        2> $PROJ_DIR/output/summaries/${SAMPLE}_summary.txt
done


## Parse summary output tot able -----
export LC_ALL=en_US.UTF-8
export PYTHONPATH=/programs/multiqc-1.10.1/lib64/python3.6/site-packages:/programs/multiqc-1.10.1/lib/python3.6/site-packages
export PATH=/programs/multiqc-1.10.1/bin:$PATH

multiqc /workdir/mbb262/nam_hybrid_rnaseq/output/summaries/


