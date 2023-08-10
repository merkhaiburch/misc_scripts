#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-12-12
# Updated... 2023-01-05
#
# Description:
# Run B73 and Mo17 through the expression alignment pipeline for Ana to test ASE
# capabilities
# ------------------------------------------------------------------------------

# Make directories
mkdir -p /workdir/mbb262/ase
mkdir /workdir/mbb262/ase/raw_reads
mkdir /workdir/mbb262/ase/trimmed_reads


## Gather data -------------------------------------------------------------------
imeta qu -d sample_title like '2021RNAHybrids_GCLBp1A05' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_0801p2E12' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_GCLTp4D7' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_GCGPp3F08' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'

# B73 are the first two
cd /workdir/mbb262/ase/raw_reads
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R087/MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R087/MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_2.fq.gz

# m017 the second two
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R088/MS21R088_CKDL210018333-2a-AK1952-AK2044_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R088/MS21R088_CKDL210018333-2a-AK1952-AK2044_HH5V7DSX2_L1_2.fq.gz

# ky21
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R021/MS21R021_CKDL210018333-2a-AK2147-7UDI223_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R021/MS21R021_CKDL210018333-2a-AK2147-7UDI223_HH5V7DSX2_L1_2.fq.gz

# additional hybrids 
#b73xmo17
imeta qu -d sample_title like '2021RNAHybrids_GCLTp4D7' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R119/MS21R119_CKDL210018333-2a-GE09-AK1954_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R119/MS21R119_CKDL210018333-2a-GE09-AK1954_HH5V7DSX2_L1_2.fq.gz

# b73 x ky21
imeta qu -d sample_title like '2021RNAHybrids_GCLTp4A6' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R108/MS21R108_CKDL210018333-2a-7UDI222-GG03_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R108/MS21R108_CKDL210018333-2a-7UDI222-GG03_HH5V7DSX2_L1_2.fq.gz

# rename for simplicity
mv MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_1.fq.gz b73_1.fq.gz
mv MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_2.fq.gz b73_2.fq.gz
mv MS21R088_CKDL210018333-2a-AK1952-AK2044_HH5V7DSX2_L1_1.fq.gz mo17_1.fq.gz
mv MS21R088_CKDL210018333-2a-AK1952-AK2044_HH5V7DSX2_L1_2.fq.gz mo17_2.fq.gz
mv MS21R119_CKDL210018333-2a-GE09-AK1954_HH5V7DSX2_L1_1.fq.gz hybrid_b73_mo17_1.fq.gz
mv MS21R119_CKDL210018333-2a-GE09-AK1954_HH5V7DSX2_L1_2.fq.gz hybrid_b73_mo17_2.fq.gz
mv MS21R021_CKDL210018333-2a-AK2147-7UDI223_HH5V7DSX2_L1_1.fq.gz ky21_1.fq.gz
mv MS21R021_CKDL210018333-2a-AK2147-7UDI223_HH5V7DSX2_L1_2.fq.gz ky21_2.fq.gz
mv MS21R108_CKDL210018333-2a-7UDI222-GG03_HH5V7DSX2_L1_1.fq.gz hybrid_b73_ky21_1.fq.gz
mv MS21R108_CKDL210018333-2a-7UDI222-GG03_HH5V7DSX2_L1_2.fq.gz hybrid_b73_ky21_2.fq.gz


## Trim reads ----------------------------------------------------------------------

# Make somewhere for the summary files to go
mkdir -p /workdir/mbb262/ase/summary_trimming

# Variables
PROJ_DIR=/workdir/mbb262/ase
RAW_READ_DIR=/workdir/mbb262/ase/raw_reads
FASTQ_TRIM=$PROJ_DIR/trimmed_reads
SUMMARY_TRIM=$PROJ_DIR/summary_trimming
N_THREADS=30

# B73
/programs/fastp-0.23.2/fastp \
  --thread $N_THREADS \
  --in1 "${RAW_READ_DIR}/b73_1.fq.gz" \
  --in2 "${RAW_READ_DIR}/b73_2.fq.gz" \
  --merge \
  --merged_out "${FASTQ_TRIM}/b73.trim.fq.gz" \
  --detect_adapter_for_pe \
  --html "${SUMMARY_TRIM}/b73.fastp.html" \
  --report_title b73

# Mo17
/programs/fastp-0.23.2/fastp \
  --thread $N_THREADS \
  --in1 "${RAW_READ_DIR}/mo17_1.fq.gz" \
  --in2 "${RAW_READ_DIR}/mo17_2.fq.gz" \
  --merge \
  --merged_out "${FASTQ_TRIM}/mo17.trim.fq.gz" \
  --detect_adapter_for_pe \
  --html "${SUMMARY_TRIM}/mo17.fastp.html" \
  --report_title mo17

# ky21
/programs/fastp-0.23.2/fastp \
  --thread $N_THREADS \
  --in1 "${RAW_READ_DIR}/ky21_1.fq.gz" \
  --in2 "${RAW_READ_DIR}/ky21_2.fq.gz" \
  --merge \
  --merged_out "${FASTQ_TRIM}/ky21.trim.fq.gz" \
  --detect_adapter_for_pe \
  --html "${SUMMARY_TRIM}/ky21.fastp.html" \
  --report_title ky21

# Hybrid - b72xmo17
/programs/fastp-0.23.2/fastp \
  --thread $N_THREADS \
  --in1 "${RAW_READ_DIR}/hybrid_b73_mo17_1.fq.gz" \
  --in2 "${RAW_READ_DIR}/hybrid_b73_mo17_2.fq.gz" \
  --merge \
  --merged_out "${FASTQ_TRIM}/hybrid_b73_mo17.trim.fq.gz" \
  --detect_adapter_for_pe \
  --html "${SUMMARY_TRIM}/hybrid_b73_mo17.fastp.html" \
  --report_title hybrid_b73_mo17

# Hybrid - b72xky21
/programs/fastp-0.23.2/fastp \
  --thread $N_THREADS \
  --in1 "${RAW_READ_DIR}/hybrid_b73_ky21_1.fq.gz" \
  --in2 "${RAW_READ_DIR}/hybrid_b73_ky21_2.fq.gz" \
  --merge \
  --merged_out "${FASTQ_TRIM}/hybrid_b73_ky21.trim.fq.gz" \
  --detect_adapter_for_pe \
  --html "${SUMMARY_TRIM}/hybrid_b73_ky21.fastp.html" \
  --report_title hybrid_b73_ky21


## Make test files for ana ----------------------------------------------------------------------------------

# Get b73 full length cDNAs: http://www.maizecdna.org/download/
# Combine them
cat maize_flcdna_1.txt maize_flcdna_group2.txt maize_flcdna_group3.txt maize_flcdna_group4.txt maize_flcdna_group5.txt > all_b73_flcdna.fasta

# Add b73 identifier
sed "/^>/s/$/_b73/" all_b73_flcdna.fasta > all_b73_flcdna_named.fasta

# Get Mo17 CDSs from reference genome
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Mo17-REFERENCE-CAU-1.0/Zm-Mo17-REFERENCE-CAU-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Mo17-REFERENCE-CAU-1.0/Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3.gz
# Only provides canonical transcripts
# Run R code here
bedtools getfasta -nameOnly -fi Zm-Mo17-REFERENCE-CAU-1.0.fa -bed mo17_cds.bed -fo mo17_cds.fasta


# Get Ky21 CDSs from reference genome
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0_Zm00031ab.1.gff3.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0_Zm00031ab.1.canonical_transcripts
# run R code here
bedtools getfasta -nameOnly -fi Zm-Ky21-REFERENCE-NAM-1.0.fa -bed ky21_cds.bed -fo ky21_cds.fasta

# Append genone names
sed "/^>/s/$/_mo17/" mo17_cds.fasta > mo17_cds_named.fasta
sed "/^>/s/$/_ky21/" ky21_cds.fasta > ky21_cds_named.fasta

# Create an artificial hybrids
cat all_b73_flcdna_named.fasta mo17_cds_named.fasta > b73_mo17_hybrid.fasta
cat all_b73_flcdna_named.fasta ky21_cds_named.fasta > b73_ky21_hybrid.fasta

# Rename files in directory
cd /workdir/mbb262/sim_reads/r_simulated/b73
mv sample_01.fasta sample_01_b73.fasta

cd /workdir/mbb262/sim_reads/r_simulated/mo17
mv sample_01.fasta sample_01_mo17.fasta

cd /workdir/mbb262/sim_reads/r_simulated/ky21
mv sample_01.fasta sample_01_ky21.fasta

cd /workdir/mbb262/sim_reads/r_simulated/hybridB73Mo17
mv sample_01.fasta sample_01_hybridB73Mo17.fasta

cd /workdir/mbb262/sim_reads/r_simulated/hybridB73Ky21
mv sample_01.fasta sample_01_hybridB73Ky21.fasta

# Copy all to a centralized directory 
cp /workdir/mbb262/sim_reads/r_simulated/b73/sample_01_b73.fasta /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads
cp /workdir/mbb262/sim_reads/r_simulated/mo17/sample_01_mo17.fasta /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads
cp /workdir/mbb262/sim_reads/r_simulated/ky21/sample_01_ky21.fasta /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads
cp /workdir/mbb262/sim_reads/r_simulated/hybridB73Mo17/sample_01_hybridB73Mo17.fasta /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads
cp /workdir/mbb262/sim_reads/r_simulated/hybridB73Ky21/sample_01_hybridB73Ky21.fasta /workdir/mbb262/sim_reads/r_simulated/copied_all_simulated_reads


## Gather and format reference transcriptomes --------------------------------------------------------------

# Get trascripts
scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/ajs692/panand_data/gene_model_annotation/fastas/Zea_mays_v5_annotatedCDS.fa /workdir/mbb262/ase

# Format the long names
sed 's/:.*//' /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS.fa > /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS_shortenedString.fa

# Subsample fasta to specific genes (in this case: canonical)
cd /workdir/mbb262/ase/references
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts
/programs/samtools-1.15.1-r/bin/samtools faidx \
    /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS_shortenedString.fa \
    -r /workdir/mbb262/ase/references/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts \
    -o /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS_shortenedString_canonical.fa

# Gather genomes
mkdir -p /workdir/mbb262/ase/references
wget -P /workdir/mbb262/ase/references https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
wget -P /workdir/mbb262/ase/references https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Mo17-REFERENCE-YAN-1.0/Zm-Mo17-REFERENCE-YAN-1.0.fa.gz
wget -P /workdir/mbb262/ase/references https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0.fa.gz

# Find matching transcript sequence in each genome by minimapping all B73 CDS sequence to each NAM genome
mkdir -p /workdir/mbb262/ase/references/nam_aligned_transcriptomes
N_THREADS=40

# b73 transcripts to b73 genome
/programs/minimap2-2.17/minimap2 -ax splice --eqx \
    -t $N_THREADS \
    -I 6G \
    /workdir/mbb262/ase/references/Zm-B73-REFERENCE-NAM-5.0.fa.gz \
    /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS_shortenedString_canonical.fa \
    > /workdir/mbb262/ase/references/nam_aligned_transcriptomes/Zm-B73-REFERENCE-NAM-5.0.sam

# b73 transcripts to mo17 genome
/programs/minimap2-2.17/minimap2 -ax splice --eqx \
    -t $N_THREADS \
    -I 6G \
    /workdir/mbb262/ase/references/Zm-Mo17-REFERENCE-YAN-1.0.fa.gz \
    /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS_shortenedString_canonical.fa \
    > /workdir/mbb262/ase/references/nam_aligned_transcriptomes/Zm-Mo17-REFERENCE-YAN-1.0.sam

# b73 transcripts to ky21 genome
/programs/minimap2-2.17/minimap2 -ax splice --eqx \
    -t $N_THREADS \
    -I 6G \
    /workdir/mbb262/ase/references/Zm-Ky21-REFERENCE-NAM-1.0.fa.gz \
    /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS_shortenedString_canonical.fa \
    > /workdir/mbb262/ase/references/nam_aligned_transcriptomes/Zm-Ky21-REFERENCE-NAM-1.0.sam



## Format reference transcriptomes ---------------------------------------------------
# (pull out sequence in each reference that aligns the best to the b73 cds seqeunce)

# Make directories
mkdir -p /workdir/mbb262/ase/references/bam_transcriptomes
mkdir -p /workdir/mbb262/ase/references/fa_transcriptomes
mkdir /workdir/mbb262/ase/references/bam_transcriptomes/primary

# Variables
SAM_TRANSCRIPTOME_DIR=/workdir/mbb262/ase/references/nam_aligned_transcriptomes
BAM_TRANSCRIPTOME_DIR=/workdir/mbb262/ase/references/bam_transcriptomes
FA_TRANSCRIPTOME_DIR=/workdir/mbb262/ase/references/fa_transcriptomes
PRIMARY_DIR=/workdir/mbb262/ase/references/bam_transcriptomes/primary

export PATH=/programs/samtools-1.15.1/bin:$PATH
export PATH=/programs/bedtools2-2.29.2/bin:$PATH

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
    /programs/samtools-1.15.1-r/bin/samtools view -e 'flag <= 16' $BAM_TRANSCRIPTOME_DIR/${SAMPLE}.bam -o $PRIMARY_DIR/primary_${SAMPLE}.bam

    # bam to bed
    bedtools bamtobed -i $PRIMARY_DIR/primary_${SAMPLE}.bam > $PRIMARY_DIR/primary_${SAMPLE}.bed
done


# Go from bed to fasta
cd /workdir/mbb262/ase/references
gunzip *.fa.gz
find /workdir/mbb262/ase/references -maxdepth 1 -name '*.fa' > $PRIMARY_DIR/refGenome.list
find $PRIMARY_DIR/ -maxdepth 1 -name '*.bed' > $PRIMARY_DIR/genomeBeds.list
cd $PRIMARY_DIR
parallel --link "bedtools getfasta -name -s -fi {1} -bed {2/} -fo {1/.}.fasta" :::: $PRIMARY_DIR/refGenome.list :::: $PRIMARY_DIR/genomeBeds.list 


# Do last bit of formatting to read names
for i in $PRIMARY_DIR/*.fasta
do
    SAMPLE=$(basename ${i} .fasta)

    # Remove where in each reference the seqeunce was pulled from
    sed 's/:.*//' $PRIMARY_DIR/${SAMPLE}.fasta > $PRIMARY_DIR/${SAMPLE}_shortNames.fasta

    # Append genome name to fasta files
    sed "/^>/s/$/_${SAMPLE}/" $PRIMARY_DIR/${SAMPLE}_shortNames.fasta > $FA_TRANSCRIPTOME_DIR/${SAMPLE}_named.fa
done


## Align back to genomes ------------------------------------------

# Create directories
mkdir -p /workdir/mbb262/ase/output/minimap_alignments

# Mapping parameters ---------------
N_THREADS=40
PROJ_DIR=/workdir/mbb262/ase
ALIGN_OUT=$PROJ_DIR/output/minimap_alignments
REF_TRANS=$PROJ_DIR/references/fa_transcriptomes
FASTQ_TRIM_MERGE_DIR=$PROJ_DIR/trimmed_reads

# Align b73 to b73
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_named.fa \
    $FASTQ_TRIM_MERGE_DIR/b73.trim.fq.gz > $ALIGN_OUT/b73_to_b73_minimap2.sam

# Align mo17 to mo17 
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    $REF_TRANS/Zm-Mo17-REFERENCE-YAN-1.0_named.fa \
    $FASTQ_TRIM_MERGE_DIR/mo17.trim.fq.gz > $ALIGN_OUT/mo17_to_mo17_minimap2.sam

# Align ky21 to ky21
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    $REF_TRANS/Zm-Ky21-REFERENCE-NAM-1.0_named.fa \
    $FASTQ_TRIM_MERGE_DIR/ky21.trim.fq.gz > $ALIGN_OUT/ky21_to_ky21_minimap2.sam

# Align hybrid b73xmo17 to b73
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_named.fa \
    $FASTQ_TRIM_MERGE_DIR/hybrid_b73_mo17.trim.fq.gz > $ALIGN_OUT/hybridB73Mo17_to_b73_minimap2.sam

# Align hybrid b73xmo17 to mo17
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    $REF_TRANS/Zm-Mo17-REFERENCE-YAN-1.0_named.fa \
    $FASTQ_TRIM_MERGE_DIR/hybrid_b73_mo17.trim.fq.gz > $ALIGN_OUT/hybridB73Mo17_to_mo17_minimap2.sam

# Align hybrid b73xky21 to b73
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_named.fa \
    $FASTQ_TRIM_MERGE_DIR/hybrid_b73_ky21.trim.fq.gz > $ALIGN_OUT/hybridB73Ky21_to_b73_minimap2.sam

# Align hybrid b73xky21 to ky21
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    $REF_TRANS/Zm-Ky21-REFERENCE-NAM-1.0_named.fa \
    $FASTQ_TRIM_MERGE_DIR/hybrid_b73_ky21.trim.fq.gz > $ALIGN_OUT/hybridB73Ky21_to_ky21_minimap2.sam


## Get mapping quality --------------------------------------------------------------

mkdir /workdir/mbb262/ase/output/minimap_alignments/minimap_stats
STAT_OUT_DIR=$PROJ_DIR/output/minimap_alignments/minimap_stats
ALIGN_OUT=$PROJ_DIR/output/minimap_alignments
N_THREADS=40

for i in $ALIGN_OUT/*.sam
do
    SAMPLE=$(basename ${i} .sam)

    echo "Getting alignment quality for: " ${SAMPLE}.sam

    /programs/samtools-1.15.1-r/bin/samtools stat \
        --threads $N_THREADS \
        ${SAMPLE}.sam > $STAT_OUT_DIR/${SAMPLE}.stat
done

# Biohpc version does not work, making my own version on conda
conda info --envs
# conda create --name py3.7 python=3.7 # unhash if conda environment is gone
conda activate py3.7
# conda install -c bioconda -c conda-forge multiqc # unhash if conda environment is gone

cd $STAT_OUT_DIR
multiqc $STAT_OUT_DIR


## Count mapped reads ------------------------------------------------------------

# Move sams into specific directories
cd $ALIGN_OUT
mkdir $ALIGN_OUT/b73
mkdir $ALIGN_OUT/mo17
mkdir $ALIGN_OUT/hybrid

mv b73_to_b73_minimap2.sam $ALIGN_OUT/b73
mv hybridB73Mo17_to_b73_minimap2.sam $ALIGN_OUT/b73
mv hybridB73Mo17_to_jointB73Mo17_minimap2.sam $ALIGN_OUT/hybrid
mv hybridB73Mo17_to_mo17_minimap2.sam $ALIGN_OUT/mo17
mv mo17_to_mo17_minimap2.sam $ALIGN_OUT/mo17

# Get transcript names
cd $REF_TRANS
grep -e ">" Zm-B73-REFERENCE-NAM-5.0_named.fa > temp.txt
sed 's/>//g' temp.txt > $ALIGN_OUT/all_b73_transcript_ids.txt

grep -e ">" Zm-Mo17-REFERENCE-YAN-1.0_named.fa > temp.txt
sed 's/>//g' temp.txt > $ALIGN_OUT/all_mo17_transcript_ids.txt

grep -e ">" b73_mo17_combined_genome.fa > temp.txt
sed 's/>//g' temp.txt > $ALIGN_OUT/all_b73_and_mo17_transcript_ids.txt

# make output dir
mkdir /workdir/mbb262/ase/output/counts

# Run Kotlin script 
# Modify in the script the name of the transcript file, the directory to the sams, and the out dir
/home/mbb262/bioinformatics/kotlinc/bin/kotlinc -script /home/mbb262/git_projects/hackathons/2022_12_compare_rna_alignment_methods/05b_test_count_noClip_sampleData.main.kts




