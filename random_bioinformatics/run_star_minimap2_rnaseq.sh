#!/usr/bin/env bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-11-12 
# Updated... 2022-04-26
# Description 
#   - Run minimap2 to align reads back to founder genome
# ---------------------------------------------------------------


# -------------------------------
#           Run STAR 
# -------------------------------


## Generate genome indexes -------

## Parameters --------------------
N_THREADS=35
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
GENOME_DIR=$PROJ_DIR/star_index/
GENOME_FA=$PROJ_DIR/references/Zm-B73-REFERENCE-NAM-5.0.fa
GTF=$PROJ_DIR/references/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3


# Generate genome indexes 
/programs/STAR/STAR \
    --runThreadN $N_THREADS \
    --runMode genomeGenerate \
    --genomeDir $GENOME_DIR \
    --genomeFastaFiles $GENOME_FA \
    --sjdbGTFfile $GTF \
    --sjdbOverhang 100


## Run STAR alignments -----------

# Parameters
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
FQ_DIR=$PROJ_DIR/reads/trimmed/
N_THREADS=35
GENOME_DIR=$PROJ_DIR/star_index/
OUT_PATH=$PROJ_DIR/output/star_alignments


# Run STAR alignment
for FQ in ${FQ_DIR}*.fq.gz
do
    FQ_BASE=$(echo "${FQ##*/}");
    /programs/STAR/STAR \
        --runThreadN $N_THREADS \
        --genomeDir $GENOME_DIR \
        --readFilesIn $FQ \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix $OUT_PATH/$FQ_BASE
done


## Organize files ----
echo "Moving files..."
mv $OUT_PATH/*.final.out $OUT_PATH/statistics/
mv $OUT_PATH/*.sam $OUT_PATH/sam/
mv $OUT_PATH/*.out.tab $OUT_PATH/splice_junctions/
rm $OUT_PATH/*.progress.out
echo "Finished!"



# ----------------------------------
#           Run minimap2
# ----------------------------------

## Gather paths --------------------
export PATH=/programs/minimap2-2.17:$PATH


## Parameters -----------------------
N_THREADS=35
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
GENOME_FA=$PROJ_DIR/references/Zm-B73-REFERENCE-NAM-5.0.fa
FASTQC_TRIM=$PROJ_DIR/reads/trimmed/*fq.gz
SAM_OUT=$PROJ_DIR/output/minimap_sam


## Run minimap2 ---------------------
minimap2 \
    -ax splice \
    $GENOME_FA \
    $FASTQC_TRIM > $SAM_OUT/${FASTQC_TRIM}_aln.sam # this outputs *fq.gz_aln.sam

