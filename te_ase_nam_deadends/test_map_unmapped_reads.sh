# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-08-11
# Updated... 2023-08-11
#
# Description:
# Select one random B73 inbred --> MS21R020
# Gather list of unmapped reads
# Count types of unpaeed reads (in R)
# Subsample those reads using seqtk subseq
# Map those reads using hitsat2
# ------------------------------------------------------------------------------

# Cerate a conda environment
conda create --name seqtk seqtk

mkdir /workdir/mbb262/temp

cp /workdir/mbb262/nam_hybrid_rnaseq/reads/trimmed/MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_trimmedQC_1.fq.gz /workdir/mbb262/temp
cp /workdir/mbb262/nam_hybrid_rnaseq/reads/trimmed/MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_trimmedQC_2.fq.gz /workdir/mbb262/temp
gunzip *.gz

head MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_trimmedQC_1.fq
head MS21R020_unmapped.lst

# subsample to unmapped reads
seqtk subseq MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_trimmedQC_1.fq MS21R020_unmapped.lst > MS21R020_1_out.fq
seqtk subseq MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_trimmedQC_2.fq MS21R020_unmapped.lst > MS21R020_2_out.fq


# Map reads using hisat2 ---------------------------------------------------------

source /programs/HISAT2/hisat2.sh


## Index parameters ----------------
mkdir /workdir/mbb262/temp/hisat2_index
mkdir /workdir/mbb262/temp/references

wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz -P /workdir/mbb262/temp/references
gunzip /workdir/mbb262/temp/references/*.gz

N_THREADS=38
PROJ_DIR=/workdir/mbb262/temp
GENOME_DIR=$PROJ_DIR/hisat2_index
GENOME_FA=$PROJ_DIR/references/Zm-B73-REFERENCE-NAM-5.0.fa


## Build index ---------------------
hisat2-build \
    $GENOME_FA \
    $GENOME_DIR


# Mapping parameters ---------------
mkdir /workdir/mbb262/temp/output/hisat2_alignments

N_THREADS=35
PROJ_DIR=/workdir/mbb262/temp
GENOME_DIR=$PROJ_DIR/hisat2_index
GENOME_FA=$PROJ_DIR/references/Zm-B73-REFERENCE-NAM-5.0.fa
OUT_PATH=$PROJ_DIR/output/hisat2_alignments


## Run hisat2 ----------------------

cd /workdir/mbb262/temp
hisat2 \
    -p $N_THREADS \
    -x $GENOME_DIR \
    -1 MS21R020_1_out.fq \
    -2 MS21R020_2_out.fq \
    -S MS21R020_output.sam --new-summary \
    2> MS21R020_summary.txt


## Parse summary output to a table -----
conda activate py3.7
multiqc /workdir/mbb262/temp/output/summaries/
conda deactivate

