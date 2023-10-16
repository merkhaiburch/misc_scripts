#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-08-10
# Updated... 2023-08-14
#
# Description:
# see spreadsheet
# ------------------------------------------------------------------------------

# # Make directories
# PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
# mkdir -p $PROJ_DIR/salmon_paired_no500_decoy_liftoffallB73
# mkdir -p $PROJ_DIR/salmon_paired_no500_decoy_liftoffallB73/index_genomes
# mkdir -p $PROJ_DIR/salmon_paired_no500_decoy_liftoffallB73/index_genomes/decoys
# mkdir -p $PROJ_DIR/salmon_paired_no500_decoy_liftoffallB73/index_genomes/combined_transcriptomes_genomes
# mkdir -p $PROJ_DIR/salmon_paired_no500_decoy_liftoffallB73/output/inbred
# mkdir -p $PROJ_DIR/salmon_paired_no500_decoy_liftoffallB73/index_genomes/liftoff

# Run liftoff ---------------------------------------------------------------------------

# Create a conda environment for this
conda create --name liftoff
conda activate liftoff
conda install -c bioconda liftoff
conda install -c bioconda agat
conda update agat

# Get genomes and gff
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz -P $PROJ_DIR/salmon_paired_no500_decoy_liftoffallB73/index_genomes/liftoff
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz -P $PROJ_DIR/salmon_paired_no500_decoy_liftoffallB73/index_genomes/liftoff/

# Unzip and bgzip
gunzip $PROJ_DIR/salmon_paired_no500_decoy_liftoffallB73/index_genomes/liftoff/Zm-B73-REFERENCE-NAM-5.0.fa.gz
bgzip --threads 20 $PROJ_DIR/salmon_paired_no500_decoy_liftoffallB73/index_genomes/liftoff/Zm-B73-REFERENCE-NAM-5.0.fa

# Run liftoff, polishing, unzipped, rearranged
cd $PROJ_DIR/salmon_paired_no500_decoy_liftoffallB73/index_genomes/liftoff2
liftoff \
    -p 45 \
    -g ./Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 \
    -o ./test_polish.gff3 \
    -polish \
    ./Zm-B73-REFERENCE-NAM-5.0_target.fa \
    ./Zm-B73-REFERENCE-NAM-5.0_reference.fa

# Run agat to get cat sequences from a gff
# https://www.biostars.org/p/9465973/
# https://github.com/NBISweden/AGAT
# https://agat.readthedocs.io/en/latest/tools/agat_sp_extract_sequences.html
agat_sp_extract_sequences.pl \
    -gff test_polish.gff3_polished \
    --fasta ./Zm-B73-REFERENCE-NAM-5.0_reference.fa \
    --mRNA \
    -o test_agat.fa

conda deactivate


# ## Make Salmon indexes for inbreds ------------------------------------------------------

# conda activate salmon

# More variables now that paths are made
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73=$PROJ_DIR/salmon_paired_no500_decoy_liftoffallB73/index_genomes/
GENOME_DIR=$PROJ_DIR/references/nam_genomes
N_THREADS=40


## Make Salmon indexes for inbreds ------------------------------------------------------

cd /workdir/mbb262/nam_hybrid_rnaseq/salmon_paired_no500_decoy_liftoffallB73/index_genomes/liftoff2

# Simplify names from fasta (remove everything after the space)
sed -i.bak -e 's/\s.*$//' test_agat.fa

# Get names for indexing
grep "^>" ./Zm-B73-REFERENCE-NAM-5.0.fa | cut -d " " -f 1 > $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/decoys/Zm-B73-REFERENCE-NAM-5.0_decoys.txt
    
# Format output to remove the carat
sed -i.bak -e 's/>//g' $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/decoys/Zm-B73-REFERENCE-NAM-5.0_decoys.txt

# Cat transcriptome and genomes
cat test_agat.fa Zm-B73-REFERENCE-NAM-5.0.fa > $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/combined_transcriptomes_genomes/Zm-B73-REFERENCE-NAM-5.0_gentrome.fa

# Index genomes
salmon index \
    --threads $N_THREADS \
    --kmerLen 31 \
    --transcripts $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/combined_transcriptomes_genomes/Zm-B73-REFERENCE-NAM-5.0_gentrome.fa \
    --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index \
    --decoys $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/decoys/Zm-B73-REFERENCE-NAM-5.0_decoys.txt


## Mapping parameters ---------------
N_THREADS=40
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
FASTQ_TRIM_DIR=$PROJ_DIR/reads/trimmed
ALIGN_OUT_LIFTOFFALLB73=$PROJ_DIR/salmon_paired_no500_decoy_liftoffallB73/output/inbred
SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73=$PROJ_DIR/salmon_paired_no500_decoy_liftoffallB73/index_genomes


## Run salmon ----------------------
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R002_CKDL210018333-2a-AK11819-AK30533_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R002_CKDL210018333-2a-AK11819-AK30533_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R002 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R003_CKDL210018333-2a-AK30535-AK30534_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R003_CKDL210018333-2a-AK30535-AK30534_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R003 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R004_CKDL210018333-2a-AK7741-AK25811_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R004_CKDL210018333-2a-AK7741-AK25811_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R004 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R005_CKDL210018333-2a-AK30537-AK30536_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R005_CKDL210018333-2a-AK30537-AK30536_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R005 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R006_CKDL210018333-2a-AK27680-AK30538_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R006_CKDL210018333-2a-AK27680-AK30538_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R006 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R007_CKDL210018333-2a-AK14231-AK30539_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R007_CKDL210018333-2a-AK14231-AK30539_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R007 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R008_CKDL210018333-2a-AK30540-AK16431_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R008_CKDL210018333-2a-AK30540-AK16431_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R008 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R009_CKDL210018333-2a-AK30542-AK30541_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R009_CKDL210018333-2a-AK30542-AK30541_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R009 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R010_CKDL210018333-2a-7UDI301-AK1847_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R010_CKDL210018333-2a-7UDI301-AK1847_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R010 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R012_CKDL210018333-2a-7UDI304-AK690_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R012_CKDL210018333-2a-7UDI304-AK690_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R012 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R013_CKDL210018333-2a-7UDI254-AK2329_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R013_CKDL210018333-2a-7UDI254-AK2329_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R013 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R015_CKDL210018333-2a-AK2044-AK1952_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R015_CKDL210018333-2a-AK2044-AK1952_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R015 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R016_CKDL210018333-2a-GH12-GF11_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R016_CKDL210018333-2a-GH12-GF11_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R016 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R017_CKDL210018333-2a-AK1868-AK1852_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R017_CKDL210018333-2a-AK1868-AK1852_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R017 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R018_CKDL210018333-2a-AK844-AK1878_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R018_CKDL210018333-2a-AK844-AK1878_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R018 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R019_CKDL210018333-2a-GG03-7UDI222_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R019_CKDL210018333-2a-GG03-7UDI222_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R019 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R020 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R021_CKDL210018333-2a-AK2147-7UDI223_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R021_CKDL210018333-2a-AK2147-7UDI223_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R021 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R022_CKDL210018333-2a-AK900-7UDI293_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R022_CKDL210018333-2a-AK900-7UDI293_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R022 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R023_CKDL210018333-2a-AK2161-AK1846_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R023_CKDL210018333-2a-AK2161-AK1846_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R023 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R024_CKDL210018333-2a-AK907-AK1864_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R024_CKDL210018333-2a-AK907-AK1864_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R024 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R025_CKDL210018333-2a-AK903-GH07_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R025_CKDL210018333-2a-AK903-GH07_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R025 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R026_CKDL210018333-2a-7UDI763-AK19481_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R026_CKDL210018333-2a-7UDI763-AK19481_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R026 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R027_CKDL210018333-2a-AK18848-AK7115_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R027_CKDL210018333-2a-AK18848-AK7115_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R027 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R028_CKDL210018333-2a-AK19520-SCI5E038_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R028_CKDL210018333-2a-AK19520-SCI5E038_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R028 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R029_CKDL210018333-2a-AK869-AK847_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R029_CKDL210018333-2a-AK869-AK847_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R029 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R034_CKDL210018333-2a-AK30548-AK30547_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R034_CKDL210018333-2a-AK30548-AK30547_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R034 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R035_CKDL210018333-2a-AK10486-AK9660_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R035_CKDL210018333-2a-AK10486-AK9660_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R035 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R036_CKDL210018333-2a-AK30550-AK30549_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R036_CKDL210018333-2a-AK30550-AK30549_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R036 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R037_CKDL210018333-2a-AK30552-AK30551_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R037_CKDL210018333-2a-AK30552-AK30551_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R037 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R038_CKDL210018333-2a-AK30553-AK15571_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R038_CKDL210018333-2a-AK30553-AK15571_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R038 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R039_CKDL210018333-2a-AK30554-AK27741_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R039_CKDL210018333-2a-AK30554-AK27741_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R039 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R040_CKDL210018333-2a-AK10434-AK30555_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R040_CKDL210018333-2a-AK10434-AK30555_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R040 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R041_CKDL210018333-2a-7UDI2717-AK23823_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R041_CKDL210018333-2a-7UDI2717-AK23823_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R041 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R042_CKDL210018333-2a-AK7701-AK30556_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R042_CKDL210018333-2a-AK7701-AK30556_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R042 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R043_CKDL210018333-2a-AK30558-AK30557_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R043_CKDL210018333-2a-AK30558-AK30557_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R043 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R044_CKDL210018333-2a-X176-AK30559_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R044_CKDL210018333-2a-X176-AK30559_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R044 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R045_CKDL210018333-2a-AK10235-AK592_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R045_CKDL210018333-2a-AK10235-AK592_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R045 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R046_CKDL210018333-2a-AK30561-AK30560_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R046_CKDL210018333-2a-AK30561-AK30560_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R046 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R047_CKDL210018333-2a-AK30562-7UDI629_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R047_CKDL210018333-2a-AK30562-7UDI629_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R047 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R048_CKDL210018333-2a-K30030-AK30563_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R048_CKDL210018333-2a-K30030-AK30563_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R048 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R049_CKDL210018333-2a-AK19627-AK30564_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R049_CKDL210018333-2a-AK19627-AK30564_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R049 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R051_CKDL210018333-2a-AK30567-AK30566_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R051_CKDL210018333-2a-AK30567-AK30566_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R051 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R052_CKDL210018333-2a-AK5231-AK4415_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R052_CKDL210018333-2a-AK5231-AK4415_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R052 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R053_CKDL210018333-2a-AK10568-AK12757_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R053_CKDL210018333-2a-AK10568-AK12757_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R053 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R054_CKDL210018333-2a-AK30569-AK30568_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R054_CKDL210018333-2a-AK30569-AK30568_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R054 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R055_CKDL210018333-2a-AK30571-AK30570_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R055_CKDL210018333-2a-AK30571-AK30570_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R055 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R056_CKDL210018333-2a-AK2209-AK15971_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R056_CKDL210018333-2a-AK2209-AK15971_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R056 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R057_CKDL210018333-2a-AK30572-AK2400_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R057_CKDL210018333-2a-AK30572-AK2400_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R057 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R058_CKDL210018333-2a-AK2486-AK30573_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R058_CKDL210018333-2a-AK2486-AK30573_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R058 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R059_CKDL210018333-2a-AK30575-AK30574_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R059_CKDL210018333-2a-AK30575-AK30574_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R059 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R060_CKDL210018333-2a-AK30576-AK7506_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R060_CKDL210018333-2a-AK30576-AK7506_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R060 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R061_CKDL210018333-2a-AK30578-AK30577_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R061_CKDL210018333-2a-AK30578-AK30577_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R061 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R062_CKDL210018333-2a-AK847-AK869_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R062_CKDL210018333-2a-AK847-AK869_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R062 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R063_CKDL210018333-2a-AK1847-7UDI301_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R063_CKDL210018333-2a-AK1847-7UDI301_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R063 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R064_CKDL210018333-2a-AK30560-AK30561_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R064_CKDL210018333-2a-AK30560-AK30561_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R064 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R065_CKDL210018333-2a-7UDI629-AK30562_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R065_CKDL210018333-2a-7UDI629-AK30562_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R065 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R066_CKDL210018333-2a-AK30563-K30030_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R066_CKDL210018333-2a-AK30563-K30030_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R066 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R067_CKDL210018333-2a-AK30564-AK19627_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R067_CKDL210018333-2a-AK30564-AK19627_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R067 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R068_CKDL210018333-2a-AK30565-AK20051_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R068_CKDL210018333-2a-AK30565-AK20051_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R068 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R069_CKDL210018333-2a-AK30566-AK30567_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R069_CKDL210018333-2a-AK30566-AK30567_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R069 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R070_CKDL210018333-2a-AK4415-AK5231_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R070_CKDL210018333-2a-AK4415-AK5231_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R070 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R071_CKDL210018333-2a-AK12757-AK10568_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R071_CKDL210018333-2a-AK12757-AK10568_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R071 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R072_CKDL210018333-2a-AK30568-AK30569_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R072_CKDL210018333-2a-AK30568-AK30569_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R072 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R073_CKDL210018333-2a-AK30570-AK30571_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R073_CKDL210018333-2a-AK30570-AK30571_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R073 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R074_CKDL210018333-2a-AK15971-AK2209_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R074_CKDL210018333-2a-AK15971-AK2209_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R074 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R075_CKDL210018333-2a-AK2400-AK30572_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R075_CKDL210018333-2a-AK2400-AK30572_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R075 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R076_CKDL210018333-2a-AK30573-AK2486_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R076_CKDL210018333-2a-AK30573-AK2486_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R076 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R078_CKDL210018333-2a-AK7506-AK30576_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R078_CKDL210018333-2a-AK7506-AK30576_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R078 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R079_CKDL210018333-2a-AK1864-AK907_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R079_CKDL210018333-2a-AK1864-AK907_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R079 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R080_CKDL210018333-2a-7UDI292-AK1334_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R080_CKDL210018333-2a-7UDI292-AK1334_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R080 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R081_CKDL210018333-2a-AK1867-7UDI237_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R081_CKDL210018333-2a-AK1867-7UDI237_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R081 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R082_CKDL210018333-2a-GD06-AK845_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R082_CKDL210018333-2a-GD06-AK845_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R082 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R083_CKDL210018333-2a-AK30577-AK30578_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R083_CKDL210018333-2a-AK30577-AK30578_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R083 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R085_CKDL210018333-2a-AK30525-AK21985_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R085_CKDL210018333-2a-AK30525-AK21985_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R085 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R087 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R088_CKDL210018333-2a-AK1952-AK2044_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R088_CKDL210018333-2a-AK1952-AK2044_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R088 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R089_CKDL210018333-2a-AK985-AK987_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R089_CKDL210018333-2a-AK985-AK987_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R089 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R090_CKDL210018333-2a-AK1872-AK898_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R090_CKDL210018333-2a-AK1872-AK898_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R090 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R092_CKDL210018333-2a-AK30529-AK30528_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R092_CKDL210018333-2a-AK30529-AK30528_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R092 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R093_CKDL210018333-2a-AK30531-AK30530_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R093_CKDL210018333-2a-AK30531-AK30530_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R093 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R094_CKDL210018333-2a-5UDI388-AK7296_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R094_CKDL210018333-2a-5UDI388-AK7296_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R094 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R095_CKDL210018333-2a-AK16359-AK30532_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R095_CKDL210018333-2a-AK16359-AK30532_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R095 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R099_CKDL210018333-2a-7UDI274-AK234_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R099_CKDL210018333-2a-7UDI274-AK234_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMON_INDEX_PAIRED_no500_DECOY_liftoffallB73/Zm-B73-REFERENCE-NAM-5.0_transcripts_index --validateMappings --writeUnmappedNames --output $ALIGN_OUT_LIFTOFFALLB73/MS21R099 


# Run multiqc -------------------------------------------

conda deactivate
conda activate py3.7
cd /workdir/mbb262/nam_hybrid_rnaseq/salmon_paired_no500_decoy_liftoffallB73/output
multiqc /workdir/mbb262/nam_hybrid_rnaseq/salmon_paired_no500_decoy_liftoffallB73/output/inbred

