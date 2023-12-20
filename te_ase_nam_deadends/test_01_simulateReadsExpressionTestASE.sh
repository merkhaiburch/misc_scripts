#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-12-12
#
# Updated... 2023-06-29
# Description:
# Simulate B73, Cml247, and Ky21 inbred and hybrid reads to test the expression 
# alignment pipeline and ASE capabilities
# ------------------------------------------------------------------------------

# Make directories
mkdir -p /workdir/mbb262/sim_reads/genomes_to_simulate_reads_from
mkdir -p /workdir/mbb262/sim_reads/r_simulated/b73
mkdir -p /workdir/mbb262/sim_reads/r_simulated/ky21
mkdir -p /workdir/mbb262/sim_reads/r_simulated/cml247
mkdir -p /workdir/mbb262/sim_reads/r_simulated/hybridB73Cml247
mkdir -p /workdir/mbb262/sim_reads/r_simulated/hybridB73Ky21
mkdir -p /workdir/mbb262/sim_reads/r_simulated/hybridCml247Ky21
mkdir /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads

cd /workdir/mbb262/sim_reads/genomes_to_simulate_reads_from

# Get B73 transcripts and best transcripts from max reelgene scores
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/te_ase_nam/top_reel_transcripts.txt ./

# Get Cml247 genome fasta and gff from reference genome, only canonical transcripts available
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML247-REFERENCE-NAM-1.0/Zm-CML247-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML247-REFERENCE-NAM-1.0/Zm-CML247-REFERENCE-NAM-1.0_Zm00023ab.1.gff3.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML247-REFERENCE-NAM-1.0/Zm-CML247-REFERENCE-NAM-1.0_Zm00023ab.1.canonical_transcripts

# Get Ky21 genome fasta and gff from reference genome, only canonical transcripts available
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0_Zm00031ab.1.gff3.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0_Zm00031ab.1.canonical_transcripts

# Unzip fa files
gunzip *.fa.gz


/programs/rstudio_server/rstudio_stop
/programs/rstudio_server/mv_dir
/programs/rstudio_server/rstudio_start 
# Run R code here



# Get fasta sequence with bedtools from canonical (Cml247, Ky21) or reelGene transcripts
bedtools getfasta -nameOnly -fi Zm-B73-REFERENCE-NAM-5.0.fa -bed b73_transcript.bed -fo b73_transcripts.fasta
bedtools getfasta -nameOnly -fi Zm-CML247-REFERENCE-NAM-1.0.fa -bed cml247_transcript.bed -fo cml247_transcripts.fasta
bedtools getfasta -nameOnly -fi Zm-Ky21-REFERENCE-NAM-1.0.fa -bed ky21_transcript.bed -fo ky21_transcripts.fasta

# Append genone names
sed "/^>/s/$/_b73/" b73_transcripts.fasta > b73_transcripts_named.fasta
sed "/^>/s/$/_cml247/" cml247_transcripts.fasta > cml247_transcripts_named.fasta
sed "/^>/s/$/_ky21/" ky21_transcripts.fasta > ky21_transcripts_named.fasta

# Create an artificial hybrids
cat b73_transcripts_named.fasta cml247_transcripts_named.fasta > b73_cml247_hybrid.fasta
cat b73_transcripts_named.fasta ky21_transcripts_named.fasta > b73_ky21_hybrid.fasta
cat cml247_transcripts_named.fasta ky21_transcripts_named.fasta > cml247_ky21_hybrid.fasta

# Simulate reads in R here ---------

# Rename files in directory (after simulating reads within R)
cd /workdir/mbb262/sim_reads/r_simulated/b73
mv sample_01.fasta sample_01_b73.fasta
mv sim_counts_matrix.rda sample_01_b73.rda

cd /workdir/mbb262/sim_reads/r_simulated/cml247
mv sample_01.fasta sample_01_cml247.fasta
mv sim_counts_matrix.rda sample_01_cml247.rda

cd /workdir/mbb262/sim_reads/r_simulated/ky21
mv sample_01.fasta sample_01_ky21.fasta
mv sim_counts_matrix.rda sample_01_ky21.rda

cd /workdir/mbb262/sim_reads/r_simulated/hybridB73Cml247
mv sample_01.fasta sample_01_hybridB73Cml247.fasta
mv sim_counts_matrix.rda sample_01_hybridB73Cml247.rda

cd /workdir/mbb262/sim_reads/r_simulated/hybridB73Ky21
mv sample_01.fasta sample_01_hybridB73Ky21.fasta
mv sim_counts_matrix.rda sample_01_hybridB73Ky21.rda

cd /workdir/mbb262/sim_reads/r_simulated/hybridCml247Ky21
mv sample_01.fasta sample_01_hybridCml247Ky21.fasta
mv sim_counts_matrix.rda sample_01_hybridCml247Ky21.rda

# Copy all fasta reads to a centralized directory, but keep a copy within each folder
cp /workdir/mbb262/sim_reads/r_simulated/b73/sample_01_b73.fasta /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads
cp /workdir/mbb262/sim_reads/r_simulated/cml247/sample_01_cml247.fasta /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads
cp /workdir/mbb262/sim_reads/r_simulated/ky21/sample_01_ky21.fasta /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads
cp /workdir/mbb262/sim_reads/r_simulated/hybridB73Cml247/sample_01_hybridB73Cml247.fasta /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads
cp /workdir/mbb262/sim_reads/r_simulated/hybridB73Ky21/sample_01_hybridB73Ky21.fasta /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads
cp /workdir/mbb262/sim_reads/r_simulated/hybridCml247Ky21/sample_01_hybridCml247Ky21.fasta /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads


## Gather and format reference transcriptomes --------------------------------------------------------------

# Make directory to store stuff
mkdir -p /workdir/mbb262/sim_reads/references
PRIMARY_DIR=/workdir/mbb262/sim_reads/references

# Get transcripts
scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/ajs692/panand_data/gene_model_annotation/fastas/Zea_mays_v5_annotatedCDS.fa $PRIMARY_DIR

# Format the long names
sed 's/:.*//' $PRIMARY_DIR/Zea_mays_v5_annotatedCDS.fa > $PRIMARY_DIR/Zea_mays_v5_annotatedCDS_shortenedString.fa

# Subsample fasta to specific genes (max reelgene)
cd $PRIMARY_DIR
/programs/samtools-1.15.1-r/bin/samtools faidx \
    $PRIMARY_DIR/Zea_mays_v5_annotatedCDS_shortenedString.fa \
    -r /workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/top_reel_transcripts.txt \
    -o $PRIMARY_DIR/Zea_mays_v5_annotatedCDS_shortenedString_reelGene.fa

# Gather genomes
wget -P $PRIMARY_DIR https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
wget -P $PRIMARY_DIR https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML247-REFERENCE-NAM-1.0/Zm-CML247-REFERENCE-NAM-1.0.fa.gz
wget -P $PRIMARY_DIR https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0.fa.gz

# Find matching transcript sequence in each genome by minimapping all B73 CDS sequence to each NAM genome
mkdir -p $PRIMARY_DIR/nam_aligned_transcriptomes
N_THREADS=40

# b73 transcripts to b73 genome
/programs/minimap2-2.17/minimap2 -ax splice --eqx \
    -t $N_THREADS \
    -I 6G \
    $PRIMARY_DIR/Zm-B73-REFERENCE-NAM-5.0.fa.gz \
    $PRIMARY_DIR/Zea_mays_v5_annotatedCDS_shortenedString_reelGene.fa \
    > $PRIMARY_DIR/nam_aligned_transcriptomes/Zm-B73-REFERENCE-NAM-5.0.sam

# b73 transcripts to cml247 genome
/programs/minimap2-2.17/minimap2 -ax splice --eqx \
    -t $N_THREADS \
    -I 6G \
    $PRIMARY_DIR/Zm-CML247-REFERENCE-NAM-1.0.fa.gz \
    $PRIMARY_DIR/Zea_mays_v5_annotatedCDS_shortenedString_reelGene.fa \
    > $PRIMARY_DIR/nam_aligned_transcriptomes/Zm-CML247-REFERENCE-NAM-1.0.sam

# b73 transcripts to ky21 genome
/programs/minimap2-2.17/minimap2 -ax splice --eqx \
    -t $N_THREADS \
    -I 6G \
    $PRIMARY_DIR/Zm-Ky21-REFERENCE-NAM-1.0.fa.gz \
    $PRIMARY_DIR/Zea_mays_v5_annotatedCDS_shortenedString_reelGene.fa \
    > $PRIMARY_DIR/nam_aligned_transcriptomes/Zm-Ky21-REFERENCE-NAM-1.0.sam


## Format reference transcriptomes ---------------------------------------------------
# (pull out sequence in each reference that aligns the best to the b73 cds seqeunce)

# Make directories
mkdir -p $PRIMARY_DIR/bam_transcriptomes
mkdir -p $PRIMARY_DIR/fa_transcriptomes
mkdir $PRIMARY_DIR/bam_transcriptomes/primary

# Variables
SAM_TRANSCRIPTOME_DIR=$PRIMARY_DIR/nam_aligned_transcriptomes
BAM_TRANSCRIPTOME_DIR=$PRIMARY_DIR/bam_transcriptomes
FA_TRANSCRIPTOME_DIR=$PRIMARY_DIR/fa_transcriptomes
PRIMARY_ALIGN_DIR=$PRIMARY_DIR/bam_transcriptomes/primary

# Loop through all transcriptomes and format them
for i in $SAM_TRANSCRIPTOME_DIR/*.sam
do
    SAMPLE=$(basename ${i} .sam)

    # echo "SAM to BAM to FA to subsampled FA: " $SAM_TRANSCRIPTOME_DIR/${SAMPLE}.sam
    echo -e "SAM to BAM to FA to FA subset: ${SAMPLE}"

    # sam to bam
    /programs/samtools-1.15.1-r/bin/samtools view -bS -@ 15 $SAM_TRANSCRIPTOME_DIR/${SAMPLE}.sam | /programs/samtools-1.15.1-r/bin/samtools sort -@ 15 -o $BAM_TRANSCRIPTOME_DIR/${SAMPLE}.bam

    # Index bam
    /programs/samtools-1.15.1-r/bin/samtools index $BAM_TRANSCRIPTOME_DIR/${SAMPLE}.bam

    # Get primary alignments from bam
    /programs/samtools-1.15.1-r/bin/samtools view -e 'flag <= 16' $BAM_TRANSCRIPTOME_DIR/${SAMPLE}.bam -o $PRIMARY_ALIGN_DIR/primary_${SAMPLE}.bam

    # bam to bed
    bedtools bamtobed -i $PRIMARY_ALIGN_DIR/primary_${SAMPLE}.bam > $PRIMARY_ALIGN_DIR/primary_${SAMPLE}.bed
done

# Go from bed to fasta
cd $PRIMARY_DIR
gunzip *0.fa.gz
find $PRIMARY_DIR -maxdepth 1 -name '*0.fa' > $PRIMARY_ALIGN_DIR/refGenome.list
find $PRIMARY_ALIGN_DIR/ -maxdepth 1 -name '*.bed' > $PRIMARY_ALIGN_DIR/genomeBeds.list
cd $FA_TRANSCRIPTOME_DIR
cd $PRIMARY_ALIGN_DIR
parallel --link "bedtools getfasta -name -s -fi {1} -bed {2/} -fo $FA_TRANSCRIPTOME_DIR/{1/.}.fasta" :::: $PRIMARY_ALIGN_DIR/refGenome.list :::: $PRIMARY_ALIGN_DIR/genomeBeds.list 

# Do last bit of formatting to read names
for i in $FA_TRANSCRIPTOME_DIR/*.fasta
do
    SAMPLE=$(basename ${i} .fasta)

    # Remove where in each reference the seqeunce was pulled from
    sed 's/:.*//' $FA_TRANSCRIPTOME_DIR/${SAMPLE}.fasta > $FA_TRANSCRIPTOME_DIR/${SAMPLE}_shortNames.fasta
done


## Align SIMULATED reads back to reference transcriptomes ------------------------------------------

# Create directories
mkdir -p /workdir/mbb262/sim_reads/output/minimap_alignments

# Mapping parameters ---------------
N_THREADS=40
PROJ_DIR=/workdir/mbb262/sim_reads
ALIGN_OUT=$PROJ_DIR/output/minimap_alignments
REF_TRANS=$PROJ_DIR/references/fa_transcriptomes
FASTQ_TRIM_MERGE_DIR=/workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads

# Align b73 to b73
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    -N 15 \
    $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta \
    $FASTQ_TRIM_MERGE_DIR/sample_01_b73.fasta > $ALIGN_OUT/b73_to_b73_minimap2.sam

# Align cml247 to cml247 
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    -N 15 \
    $REF_TRANS/Zm-CML247-REFERENCE-NAM-1.0_shortNames.fasta \
    $FASTQ_TRIM_MERGE_DIR/sample_01_cml247.fasta > $ALIGN_OUT/cml247_to_cml247_minimap2.sam

# Align ky21 to ky21
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    -N 15 \
    $REF_TRANS/Zm-Ky21-REFERENCE-NAM-1.0_shortNames.fasta \
    $FASTQ_TRIM_MERGE_DIR/sample_01_ky21.fasta > $ALIGN_OUT/ky21_to_ky21_minimap2.sam

# Align hybrid b73xcml247 to b73
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    -N 15 \
    $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta \
    $FASTQ_TRIM_MERGE_DIR/sample_01_hybridB73Cml247.fasta > $ALIGN_OUT/hybridB73Cml247_to_b73_minimap2.sam

# Align hybrid b73xcml247 to cml247
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    -N 15 \
    $REF_TRANS/Zm-CML247-REFERENCE-NAM-1.0_shortNames.fasta \
    $FASTQ_TRIM_MERGE_DIR/sample_01_hybridB73Cml247.fasta> $ALIGN_OUT/hybridB73Cml247_to_cml247_minimap2.sam

# Align hybrid b73xky21 to b73
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    -N 15 \
    $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta \
    $FASTQ_TRIM_MERGE_DIR/sample_01_hybridB73Ky21.fasta > $ALIGN_OUT/hybridB73Ky21_to_b73_minimap2.sam

# Align hybrid b73xky21 to ky21
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    -N 15 \
    $REF_TRANS/Zm-Ky21-REFERENCE-NAM-1.0_shortNames.fasta \
    $FASTQ_TRIM_MERGE_DIR/sample_01_hybridB73Ky21.fasta > $ALIGN_OUT/hybridB73Ky21_to_ky21_minimap2.sam

# Align hybrid cml247xky21 to cml247
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    -N 15 \
    $REF_TRANS/Zm-CML247-REFERENCE-NAM-1.0_shortNames.fasta \
    $FASTQ_TRIM_MERGE_DIR/sample_01_hybridCml247Ky21.fasta> $ALIGN_OUT/hybridCml247Ky21_to_cml247_minimap2.sam

# Align hybrid cml247xky21 to ky21
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    -N 15 \
    $REF_TRANS/Zm-Ky21-REFERENCE-NAM-1.0_shortNames.fasta \
    $FASTQ_TRIM_MERGE_DIR/sample_01_hybridCml247Ky21.fasta > $ALIGN_OUT/hybridCml247Ky21_to_ky21_minimap2.sam


## Get mapping quality --------------------------------------------------------------

mkdir /workdir/mbb262/sim_reads/output/minimap_alignments/minimap_stats
STAT_OUT_DIR=$PROJ_DIR/output/minimap_alignments/minimap_stats
ALIGN_OUT=$PROJ_DIR/output/minimap_alignments
N_THREADS=40

cd $ALIGN_OUT
for i in $ALIGN_OUT/*.sam
do
    SAMPLE=$(basename ${i} .sam)

    echo "Getting alignment quality for: " ${SAMPLE}.sam

    /programs/samtools-1.15.1-r/bin/samtools stat \
        --threads $N_THREADS \
        ${SAMPLE}.sam > $STAT_OUT_DIR/${SAMPLE}.stat
done

# Biohpc version does not work, making my own version on conda
# conda info --envs
# conda create --name py3.7 python=3.7 # unhash if conda environment is gone
conda activate py3.7
# conda install -c bioconda -c conda-forge multiqc # unhash if conda environment is gone

cd $STAT_OUT_DIR
multiqc $STAT_OUT_DIR


## Count mapped reads ------------------------------------------------------------

## Set variables ---------------------

# minimum alignment length proportion
MIN_ALIGN=0.9

# minimum mapQ score
MIN_MAPQ=48

# "relaxed" maximum edit distance proportion
# used to find subpar matches to both parents
LAX_NM=0.1

# "strict" maximum edit distance proportion
# only reads that pass this threshold on at least one parent are counted
STRICT_NM=0.02

# Inbred maximum edit distance proportion
INBRED_NM=0.02


## Files ------------------------------

mkdir /workdir/mbb262/sim_reads/output/read_counts
PATH_TO_SAMS=/workdir/mbb262/sim_reads/output/minimap_alignments
PATH_TO_OUT_COUNTS=/workdir/mbb262/sim_reads/output/read_counts

SAM_2="hybridB73Ky21_to_b73_minimap2.sam"
SAM_3="hybridB73Ky21_to_ky21_minimap2.sam"
SAM_4="hybridB73Cml247_to_b73_minimap2.sam"
SAM_5="hybridB73Cml247_to_cml247_minimap2.sam"
SAM_6="hybridCml247Ky21_to_cml247_minimap2.sam"
SAM_7="hybridCml247Ky21_to_ky21_minimap2.sam"

# suffixes are the reference name added to the end of the transcript name
# for each transcript/contig
SUFFIX_2="_b73"
SUFFIX_3="_ky21"
SUFFIX_4="_b73"
SUFFIX_5="_cml247"
SUFFIX_6="_cml247"
SUFFIX_7="_ky21"


# file containing a list of transcript id's to count
# e.g. Zm00001eb233020_T002
# one per line of the file
TRANSCRIPT_IDS="/workdir/mbb262/sim_reads/r_simulated/simulated_transcript_ids.txt"

# output file
OUTPUT_1="inbred_counts_minimap2.txt"
OUTPUT_2="hybridB73Ky21_minimap2.txt"
OUTPUT_3="hybridB73Cml247_minimap2.txt"
OUTPUT_4="hybridCml247Ky21_minimap2.txt"

## run the scripts -----------------------

# count_rnaseq_reads_hybrid.main.kts --> ASE counter (Ana)
# 04b_count_rnaseq_reads_minimap2.main.kts --> inbred counter (Zack)

# B73 to B73, Ky21 to Ky21, and Cml247 to Cml247 (does everything in a directory)
/home/mbb262/bioinformatics/kotlinc/bin/kotlinc \
    -script /home/mbb262/git_projects/te_ase_nam/src/expression_04b_countRNAseqReadsMinimap2Script.main.kts -- \
    -i $PATH_TO_SAMS \
    -t $TRANSCRIPT_IDS \
    -o $PATH_TO_OUT_COUNTS/$OUTPUT_1 \
    -a $MIN_ALIGN \
    -m $INBRED_NM \
    -n $MIN_MAPQ

# hybrid B73 x Ky21
/home/mbb262/bioinformatics/kotlinc/bin/kotlinc \
    -script /home/mbb262/git_projects/te_ase_nam/src/count_rnaseq_reads_hybrid.main.kts -- \
    -i $PATH_TO_SAMS/$SAM_2 \
    -j $PATH_TO_SAMS/$SAM_3 \
    -p $SUFFIX_2 \
    -q $SUFFIX_3 \
    -t $TRANSCRIPT_IDS \
    -o $PATH_TO_OUT_COUNTS/$OUTPUT_2 \
    -a $MIN_ALIGN \
    -l $LAX_NM \
    -m $STRICT_NM \
    -n $MIN_MAPQ

# hybrid B73 x Cml247
/home/mbb262/bioinformatics/kotlinc/bin/kotlinc \
    -script /home/mbb262/git_projects/te_ase_nam/src/count_rnaseq_reads_hybrid.main.kts -- \
    -i $PATH_TO_SAMS/$SAM_4 \
    -j $PATH_TO_SAMS/$SAM_5 \
    -p $SUFFIX_4 \
    -q $SUFFIX_5 \
    -t $TRANSCRIPT_IDS \
    -o $PATH_TO_OUT_COUNTS/$OUTPUT_3 \
    -a $MIN_ALIGN \
    -l $LAX_NM \
    -m $STRICT_NM \
    -n $MIN_MAPQ

# hybrid Cml247 x Ky21
/home/mbb262/bioinformatics/kotlinc/bin/kotlinc \
    -script /home/mbb262/git_projects/te_ase_nam/src/count_rnaseq_reads_hybrid.main.kts -- \
    -i $PATH_TO_SAMS/$SAM_6 \
    -j $PATH_TO_SAMS/$SAM_7 \
    -p $SUFFIX_6 \
    -q $SUFFIX_7 \
    -t $TRANSCRIPT_IDS \
    -o $PATH_TO_OUT_COUNTS/$OUTPUT_4 \
    -a $MIN_ALIGN \
    -l $LAX_NM \
    -m $STRICT_NM \
    -n $MIN_MAPQ


# Download maize pan-genome table
cd /workdir/mbb262/sim_reads/r_simulated
wget https://de.cyverse.org/anon-files//iplant/home/shared/NAM/NAM_genome_and_annotation_Jan2021_release/SUPPLEMENTAL_DATA/pangene-files/pan_gene_matrix_v3_cyverse.csv


# ---------------------------------------------------------------------------
## Try using salmon and seesaw to count reads and ASE reads
# ---------------------------------------------------------------------------

mkdir /workdir/mbb262/sim_reads/salmon
mkdir /workdir/mbb262/sim_reads/salmon/b73

REF_TRANS=$PROJ_DIR/references/fa_transcriptomes
$REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta

FASTQ_TRIM_MERGE_DIR=/workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads

$FASTQ_TRIM_MERGE_DIR/sample_01_b73.fasta

# create a conda environment
# conda create --name salmon
conda activate salmon
conda install -c bioconda salmon

# Salmon for inbred b73 -----------------------------

# Make a decoy 
cd /workdir/mbb262/sim_reads/salmon/b73
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz

# Get names for indexing
# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
grep "^>" <(gunzip -c Zm-B73-REFERENCE-NAM-5.0.fa.gz) | cut -d " " -f 1 > b73_decoys.txt
sed -i.bak -e 's/>//g' b73_decoys.txt

# Copy over transcripts from other directory
cp /workdir/mbb262/sim_reads/references/fa_transcriptomes/Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta ./

# cat transcriptome and genome
cat Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta <(bgzip --decompress --stdout Zm-B73-REFERENCE-NAM-5.0.fa.gz) > b73_gentrome.fa

# index
salmon index \
    --threads 40 \
    --kmerLen 31 \
    --transcripts b73_gentrome.fa \
    --index b73_transcripts_index \
    --decoys b73_decoys.txt

# quantify
salmon quant \
    --threads 40 \
    --libType A \
    -r /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads/sample_01_b73.fasta \
    --index b73_transcripts_index \
    --validateMappings \
    --output b73_transcripts_quant


# Salmon for inbred cml247 -----------------------------

# Make a decoy 
mkdir /workdir/mbb262/sim_reads/salmon/cml247
cd /workdir/mbb262/sim_reads/salmon/cml247
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML247-REFERENCE-NAM-1.0/Zm-CML247-REFERENCE-NAM-1.0.fa.gz

# Get names for indexing
# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
grep "^>" <(gunzip -c Zm-CML247-REFERENCE-NAM-1.0.fa.gz) | cut -d " " -f 1 > cml247_decoys.txt
sed -i.bak -e 's/>//g' cml247_decoys.txt

# Copy over transcripts from other directory
cp /workdir/mbb262/sim_reads/references/fa_transcriptomes/Zm-CML237-REFERENCE-NAM-5.0_shortNames.fasta ./

# cat transcriptome and genome
cat Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta <(bgzip --decompress --stdout Zm-CML247-REFERENCE-NAM-1.0.fa.gz) > cml247_gentrome.fa

# index
salmon index \
    --threads 40 \
    --kmerLen 31 \
    --transcripts cml247_gentrome.fa \
    --index cml247_transcripts_index \
    --decoys cml247_decoys.txt

# quantify
salmon quant \
    --threads 40 \
    --libType A \
    -r /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads/sample_01_cml247.fasta \
    --index cml247_transcripts_index \
    --validateMappings \
    --output cml247_transcripts_quant


# Salmon for inbred ky21 -----------------------------

# Make a decoy 
mkdir /workdir/mbb262/sim_reads/salmon/ky21
cd /workdir/mbb262/sim_reads/salmon/ky21
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0.fa.gz

# Get names for indexing
# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
grep "^>" <(gunzip -c Zm-Ky21-REFERENCE-NAM-1.0.fa.gz) | cut -d " " -f 1 > ky21_decoys.txt
sed -i.bak -e 's/>//g' ky21_decoys.txt

# Copy over transcripts from other directory
cp /workdir/mbb262/sim_reads/references/fa_transcriptomes/Zm-Ky21-REFERENCE-NAM-1.0_shortNames.fasta ./

# cat transcriptome and genome
cat Zm-Ky21-REFERENCE-NAM-1.0_shortNames.fasta <(bgzip --decompress --stdout Zm-Ky21-REFERENCE-NAM-1.0.fa.gz) > ky21_gentrome.fa

# index
salmon index \
    --threads 40 \
    --kmerLen 31 \
    --transcripts ky21_gentrome.fa \
    --index ky21_transcripts_index \
    --decoys ky21_decoys.txt

# quantify
salmon quant \
    --threads 40 \
    --libType A \
    -r /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads/sample_01_ky21.fasta \
    --index ky21_transcripts_index \
    --validateMappings \
    --output ky21_transcripts_quant


# Salmon for hybrids -----------------------------

# Make a decoy 
mkdir /workdir/mbb262/sim_reads/salmon/hybridb73ky21
cd /workdir/mbb262/sim_reads/salmon/hybridb73ky21
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0.fa.gz
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz

# Add a identifier with the genome name
sed "/^>/s/$/_ky21/" <(gunzip -c Zm-Ky21-REFERENCE-NAM-1.0.fa.gz) > ky21_genome_labeled.fa
sed "/^>/s/$/_b73/" <(gunzip -c Zm-B73-REFERENCE-NAM-5.0.fa.gz) > b73_genome_labeled.fa

# Get genome chromosome and scaffold names for indexing
# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
grep "^>" ky21_genome_labeled.fa | cut -d " " -f 1 > ky21_decoys.txt
grep "^>" b73_genome_labeled.fa | cut -d " " -f 1 > b73_decoys.txt

# Combine decoys and remove greater than symbols
cat ky21_decoys.txt b73_decoys.txt > b73_ky21_joint_decoys.txt
sed -i.bak -e 's/>//g' b73_ky21_joint_decoys.txt

# Merge transcritomes over genome from other directory after appending genome name
sed "/^>/s/$/_b73/" /workdir/mbb262/sim_reads/references/fa_transcriptomes/Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta > Zm-B73-REFERENCE-NAM-5.0_shortNames_genome_labeled.fa
sed "/^>/s/$/_ky21/" /workdir/mbb262/sim_reads/references/fa_transcriptomes/Zm-Ky21-REFERENCE-NAM-1.0_shortNames.fasta > Zm-Ky21-REFERENCE-NAM-1.0_shortNames_genome_labeled.fa

cat Zm-Ky21-REFERENCE-NAM-1.0_shortNames_genome_labeled.fa \
    Zm-B73-REFERENCE-NAM-5.0_shortNames_genome_labeled.fa \
    ky21_genome_labeled.fa \
    b73_genome_labeled.fa > b73ky21_gentrome.fa

# generate diploid txome with g2gtools:
# http://churchill-lab.github.io/g2gtools/
salmon index \
    --threads 40 \
    --kmerLen 31 \
    --transcripts b73ky21_gentrome.fa \
    --index diploid_txome_transcripts_index \
    --decoys b73_ky21_joint_decoys.txt \
    --keepDuplicates

salmon quant \
    -i diploid_txome_transcripts_index \
    --libType A \
    -r /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads/sample_01_hybridB73Ky21.fasta \
    --threads 40 \
    --numBootstraps 30 \
    -o hybrid_transcripts_quant




salmon quant \
    --threads $N_THREADS \
    --libType A \
    -k 31 \
    -1 $FASTQ_TRIM_MERGE_DIR/${sample_1} \
    -2 \$FASTQ_TRIM_MERGE_DIR/${sample_2} \
    --index \$SALMON_INDEX/$inbred1 \
    --validateMappings \
    --writeUnmappedNames \
    --output \$ALIGN_OUT/$outid













