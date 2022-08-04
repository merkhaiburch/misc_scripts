# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-04-27 
# Updated... 2022-06-28
#
# Description:
# - Test script to run hisat to align reads back to their super 
#   founder genome using hisat 2
# ---------------------------------------------------------------

# Gather genomes -------------------------------------------------------------------------------

# make directories
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/references
cd /workdir/mbb262/nam_hybrid_rnaseq/references

# Gather b73 
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz

# Gather Ky21
wget https://download.maizegdb.org/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0.fa.gz
wget https://download.maizegdb.org/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0_Zm00031ab.1.gff3.gz

# unzip files
gunzip /workdir/mbb262/nam_hybrid_rnaseq/references/*.gz

# Append the genome's name to each chromosome/scaffold name to keep track of which genome is which
sed '/^>/s/$/_b73/' Zm-B73-REFERENCE-NAM-5.0.fa > Zm-B73-REFERENCE-NAM-5.0_named.fa
sed '/^>/s/$/_ky21/' Zm-Ky21-REFERENCE-NAM-1.0.fa > Zm-Ky21-REFERENCE-NAM-1.0_named.fa


# Combine genomes + gffs ---------------------------------------------------------------------------

cat Zm-B73-REFERENCE-NAM-5.0_named.fa Zm-Ky21-REFERENCE-NAM-1.0_named.fa > b73_ky21.fa

# convert gff3 to gtf --> b73
export PATH=/programs/cufflinks-2.2.1.Linux_x86_64:$PATH
gffread \
    /workdir/mbb262/nam_hybrid_rnaseq/references/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 \
    -T -o \
    /workdir/mbb262/nam_hybrid_rnaseq/references/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gtf

# convert gff3 to gtf --> ky21
export PATH=/programs/cufflinks-2.2.1.Linux_x86_64:$PATH
gffread \
    /workdir/mbb262/nam_hybrid_rnaseq/references/Zm-Ky21-REFERENCE-NAM-1.0_Zm00031ab.1.gff3 \
    -T -o \
    /workdir/mbb262/nam_hybrid_rnaseq/references/Zm-Ky21-REFERENCE-NAM-1.0_Zm00031ab.1.gtf


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


# Run TrimGalore ---------------------------------------------------------------------------------

# Make directories
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/reads/raw
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/reads/trimmed
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/output/fastqc

## Parameters ----
N_THREADS=38
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
FASTQC_RAW_DIR=$PROJ_DIR/reads/raw 
FASTQC_OUT=$PROJ_DIR/output/fastqc/
FASTQC_TRIM=$PROJ_DIR/reads/trimmed/


## Trim Galore! for paired end reads -------------------------------------------------------------
for i in $FASTQC_RAW_DIR/*_1.fq.gz
do
    SAMPLE=$(basename ${i} _1.fq.gz)

    echo "Trimming: " ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz

    trim_galore \
    --cores $N_THREADS \
    --paired $FASTQC_RAW_DIR/${SAMPLE}_1.fq.gz $FASTQC_RAW_DIR/${SAMPLE}_2.fq.gz \
    --output_dir $FASTQC_TRIM \
    --quality 20 \
    --fastqc \
    --fastqc_args "-o $FASTQC_OUT" \
    --basename $SAMPLE
done

# move summary outputs to a differnt directory
mkdir /workdir/mbb262/nam_hybrid_rnaseq/reads/trimming_reports
mv /workdir/mbb262/nam_hybrid_rnaseq/reads/trimmed/*_trimming_report.txt /workdir/mbb262/nam_hybrid_rnaseq/reads/trimming_reports


# Run Hisat2 on super genome -----------------------------------------------------------------------

## Gather paths --------------------
source /programs/HISAT2/hisat2.sh

## Index parameters ----------------
mkdir /workdir/mbb262/nam_hybrid_rnaseq/hisat2_index
N_THREADS=38
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
GENOME_DIR=$PROJ_DIR/hisat2_index
GENOME_FA=$PROJ_DIR/references/b73_ky21.fa


## Build index ---------------------
cd $GENOME_DIR
hisat2-build \
    $GENOME_FA \
    $GENOME_DIR

# b73 only
hisat2-build \
    /workdir/mbb262/nam_hybrid_rnaseq/references/Zm-B73-REFERENCE-NAM-5.0.fa \
    /workdir/mbb262/nam_hybrid_rnaseq/hisat_index_b73
mv /workdir/mbb262/nam_hybrid_rnaseq/hisat_index_b73.* /workdir/mbb262/nam_hybrid_rnaseq/hisat_index_b73

# ky21
hisat2-build \
    /workdir/mbb262/nam_hybrid_rnaseq/references/Zm-Ky21-REFERENCE-NAM-1.0.fa \
    /workdir/mbb262/nam_hybrid_rnaseq/hisat_index_ky21
mv /workdir/mbb262/nam_hybrid_rnaseq/hisat_index_ky21.* /workdir/mbb262/nam_hybrid_rnaseq/hisat_index_ky21


# Mapping parameters ---------------
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/output/hisat2_alignments
mkdir /workdir/mbb262/nam_hybrid_rnaseq/output/summaries

source /programs/HISAT2/hisat2.sh
N_THREADS=25
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
FASTQC_TRIM_DIR=$PROJ_DIR/reads/trimmed
GENOME_DIR=$PROJ_DIR/hisat2_index
OUT_PATH=$PROJ_DIR/output/hisat2_alignments


## Run hisat2 -------------------------------------------------------------------------------------

# Map all reads to b73xky21 --> k20
for i in $FASTQC_TRIM_DIR/*_1.fq.gz
do
    SAMPLE=$(basename ${i} _1.fq.gz)

    echo "Aligning: " ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz

    hisat2 \
        -p $N_THREADS \
        -x $GENOME_DIR/hisat2_index \
        -k 20 \
        --no-softclip \
        -1 $FASTQC_TRIM_DIR/${SAMPLE}_1.fq.gz \
        -2 $FASTQC_TRIM_DIR/${SAMPLE}_2.fq.gz \
        -S $OUT_PATH/${SAMPLE}_b732ky21_k20_output.sam --new-summary \
        2> $PROJ_DIR/output/summaries/${SAMPLE}_b732ky21_k20_summary.txt
done
