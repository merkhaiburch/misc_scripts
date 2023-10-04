# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-09-21 
# Updated... 2023-08-09
#
# Description:
# Align all B73 transcripts with minimap2 to each NAM founder genome, fish out
# that sequence to be used by salmon in the next script
# ------------------------------------------------------------------------------

## Make directories ------------------------------------------------------------
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq

mkdir -p $PROJ_DIR/references
mkdir -p $PROJ_DIR/references/nam_genomes
mkdir -p $PROJ_DIR/references/bam_transcriptomes
mkdir -p $PROJ_DIR/references/fa_transcriptomes
mkdir -p $PROJ_DIR/references/nam_aligned_transcriptomes
mkdir -p $PROJ_DIR/references/bam_transcriptomes/primary

# Set variables
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
NAM_GENOME_DIR=$PROJ_DIR/references/nam_genomes
SAM_TRANSCRIPTOME_DIR=$PROJ_DIR/references/nam_aligned_transcriptomes
BAM_TRANSCRIPTOME_DIR=$PROJ_DIR/references/bam_transcriptomes
FA_TRANSCRIPTOME_DIR=$PROJ_DIR/references/fa_transcriptomes
PRIMARY_ALIGN_DIR=$PROJ_DIR/references/bam_transcriptomes/primary

N_THREADS=40


## Gather transcripts ----------------------------------------------------------

# Get list of transcripts with the max reelGene score (calculated in *top_reel_transcripts.R)
scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/te_ase_nam/top_reel_transcripts.txt $NAM_GENOME_DIR

# Download B73 v5 genome and gff
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz -P $NAM_GENOME_DIR
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz -P $NAM_GENOME_DIR
gunzip $NAM_GENOME_DIR/Zm-B73-REFERENCE-NAM-5.0.fa.gz
gunzip $NAM_GENOME_DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz

# Extract coordinates of all transcripts from v5 gff 
# Add on 500 bp to the starts and stops of each transcript to capture more UTR
# The UTRs are supposed to be included in the mRNAs but Ed said to do this because of reelGene
cd $NAM_GENOME_DIR
grep $'\tmRNA\t' Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 | awk -F '\t' -v OFS='\t' '{split($9, attributes, ";"); transcript_id = ""; for (i in attributes) {if (attributes[i] ~ /^ID=/) {transcript_id = substr(attributes[i], 4); break;}} print $1, $4 - 501, $5 + 500, transcript_id;}' > mRNA_bedfile.bed
# grep $'\tmRNA\t' Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 | awk -F '\t' -v OFS='\t' '{split($9, attributes, ";"); transcript_id = ""; for (i in attributes) {if (attributes[i] ~ /^ID=/) {transcript_id = substr(attributes[i], 4); break;}} print $1, $4 - 1, $5, transcript_id;}' > mRNA_bedfile.bed

# Remove lines with negative start sites
sed -i -r '/-/d' $NAM_GENOME_DIR/mRNA_bedfile.bed

# Get fasta sequences of all transcripts from gff 
# skips chromosomes and scaffolds where end is beyond scaffold length
bedtools getfasta -nameOnly \
    -fi $NAM_GENOME_DIR/Zm-B73-REFERENCE-NAM-5.0.fa \
    -bed $NAM_GENOME_DIR/mRNA_bedfile.bed \
    -fo $NAM_GENOME_DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.mRNA.fa

# Subsample fa file to specific genes (max reelgene)
/programs/samtools-1.15.1-r/bin/samtools faidx \
    $NAM_GENOME_DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.mRNA.fa \
    -r $NAM_GENOME_DIR/top_reel_transcripts.txt \
    -o $NAM_GENOME_DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.mRNA_reelGene.fasta

# Remove intermediate files
rm $NAM_GENOME_DIR/mRNA_bedfile.bed
rm $NAM_GENOME_DIR/top_reel_transcripts.txt
rm $NAM_GENOME_DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3
rm $NAM_GENOME_DIR/Zm-B73-REFERENCE-NAM-5.0.f*
rm $NAM_GENOME_DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.mRNA.fa
rm $NAM_GENOME_DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.mRNA.fa.fai


## Gather NAM genomes -----------------------------------------------------------

# Change directory to where all the genomes should be stored
cd $NAM_GENOME_DIR

# PHZ51 is on blfs1, change suffix
scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/jav246/cmlAssemblies/assemblyCollection/5_Buckler-PHZ51_mecatErrorCorrected.contigs.fasta.gz ./
mv 5_Buckler-PHZ51_mecatErrorCorrected.contigs.fasta.gz 5_Buckler-PHZ51_mecatErrorCorrected.contigs.fa.gz

# PHB47 is available through UGA, download it manually and upload it to cbsu
# https://acsess.onlinelibrary.wiley.com/doi/10.1002/tpg2.20114
# http://maize.uga.edu/maize_stiff_stalk_download.shtml
# PHB47 v1 genome assembly
# Change suffix after downloading
wget http://maize.uga.edu/data/PHB47_v1.asm.fa.gz

# Download remaining NAM genomes
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B97-REFERENCE-NAM-1.0/Zm-B97-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML322-REFERENCE-NAM-1.0/Zm-CML322-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML333-REFERENCE-NAM-1.0/Zm-CML333-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML247-REFERENCE-NAM-1.0/Zm-CML247-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML228-REFERENCE-NAM-1.0/Zm-CML228-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML277-REFERENCE-NAM-1.0/Zm-CML277-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML103-REFERENCE-NAM-1.0/Zm-CML103-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML69-REFERENCE-NAM-1.0/Zm-CML69-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML52-REFERENCE-NAM-1.0/Zm-CML52-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-HP301-REFERENCE-NAM-1.0/Zm-HP301-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Il14H-REFERENCE-NAM-1.0/Zm-Il14H-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Ki3-REFERENCE-NAM-1.0/Zm-Ki3-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Ki11-REFERENCE-NAM-1.0/Zm-Ki11-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-M37W-REFERENCE-NAM-1.0/Zm-M37W-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-M162W-REFERENCE-NAM-1.0/Zm-M162W-REFERENCE-NAM-1.0.fa.gz
wget https://download.maizegdb.org/Zm-Mo17-REFERENCE-CAU-2.0/Zm-Mo17-REFERENCE-CAU-2.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Mo18W-REFERENCE-NAM-1.0/Zm-Mo18W-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Ms71-REFERENCE-NAM-1.0/Zm-Ms71-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-NC350-REFERENCE-NAM-1.0/Zm-NC350-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-NC358-REFERENCE-NAM-1.0/Zm-NC358-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Oh7B-REFERENCE-NAM-1.0/Zm-Oh7B-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Oh43-REFERENCE-NAM-1.0/Zm-Oh43-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-P39-REFERENCE-NAM-1.0/Zm-P39-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Tx303-REFERENCE-NAM-1.0/Zm-Tx303-REFERENCE-NAM-1.0.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Tzi8-REFERENCE-NAM-1.0/Zm-Tzi8-REFERENCE-NAM-1.0.fa.gz


## Align B73  transcripts to each of the NAM founders -----------------------------

# Create loop for all genomes
for i in $NAM_GENOME_DIR/*.fa.gz
do
    SAMPLE=$(basename ${i} .fa.gz)

    echo "Aligning: " ${SAMPLE}.fa.gz

    /programs/minimap2-2.17/minimap2 -ax splice --eqx \
        -t $N_THREADS \
        -I 6G \
        $NAM_GENOME_DIR/${SAMPLE}.fa.gz \
        $NAM_GENOME_DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.mRNA_reelGene.fasta \
        > $SAM_TRANSCRIPTOME_DIR/${SAMPLE}_axSplice_eqx.sam
done


## Save to blfs1 ------------------------------------------------------------------

scp $PROJ_DIR/references/nam_aligned_transcriptomes/*_axSplice_eqx.sam mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/te_ase_nam/nam_transcriptomes_primary

# Get sam files if previosuly generated
# cd $PROJ_DIR/references/nam_aligned_transcriptomes
# scp -r mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/te_ase_nam/nam_transcriptomes_primary/* $PROJ_DIR/references/nam_aligned_transcriptomes/


## Format SAM files to bed files ---------------------------------------------------

# Variables for converting between file types

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
    bedtools bamtobed -i $PRIMARY_ALIGN_DIR/primary_${SAMPLE}.bam > $PRIMARY_ALIGN_DIR/${SAMPLE}_primary.bed

done


## From bed files get fasta --------------------------------------------------------

# Move to directory and unzip files
cd $NAM_GENOME_DIR
find $NAM_GENOME_DIR/*.gz -type f > $NAM_GENOME_DIR/unzippable_files.list
/programs/parallel/bin/parallel -j 11 "gunzip {}" :::: $NAM_GENOME_DIR/unzippable_files.list

# Find all unzipped reference genomes & sort alphabeticially
find $NAM_GENOME_DIR -maxdepth 1 -name '*.fa' > $PRIMARY_ALIGN_DIR/refGenome_unsorted.list
sort $PRIMARY_ALIGN_DIR/refGenome_unsorted.list > $PRIMARY_ALIGN_DIR/refGenome.list

# Find all bed files indicating fished out sequence coordinates & sort alphabeticially
find $PRIMARY_ALIGN_DIR/ -maxdepth 1 -name '*.bed' > $PRIMARY_ALIGN_DIR/genomeBeds_unsorted.list
sort $PRIMARY_ALIGN_DIR/genomeBeds_unsorted.list > $PRIMARY_ALIGN_DIR/genomeBeds.list

# Check that both files match up (in terms of genomes listed in same order)
head $PRIMARY_ALIGN_DIR/refGenome.list
head $PRIMARY_ALIGN_DIR/genomeBeds.list

# Move to correct directory and run bedtools to getfasta
cd $PRIMARY_ALIGN_DIR
parallel --link "bedtools getfasta -name -s -fi {1} -bed {2/} -fo {1/.}.fasta" :::: $PRIMARY_ALIGN_DIR/refGenome.list :::: $PRIMARY_ALIGN_DIR/genomeBeds.list 

# Do last bit of formatting to read names to remove all extra alignment names and positions
# Just keeps gene and transcript name
for i in $PRIMARY_ALIGN_DIR/*.fasta
do
    SAMPLE=$(basename ${i} .fasta)

    # Remove where in each reference the seqeunce was pulled from
    sed 's/:.*//' $PRIMARY_ALIGN_DIR/${SAMPLE}.fasta > $FA_TRANSCRIPTOME_DIR/${SAMPLE}_shortNames.fa
done



## Get all names of files for mapping script ------------------------------------------

ll $FA_TRANSCRIPTOME_DIR > /workdir/mbb262/genome_names.txt


## Get maping quality for aligned transcripts -----------------------------------------
cd $SAM_TRANSCRIPTOME_DIR
mkdir $SAM_TRANSCRIPTOME_DIR/stats

for i in $SAM_TRANSCRIPTOME_DIR/*.sam
do
    SAMPLE=$(basename ${i} .sam)

    echo "Getting alignment quality for: " ${SAMPLE}.sam

    /programs/samtools-1.15.1-r/bin/samtools stat \
        --threads $N_THREADS \
        $SAM_TRANSCRIPTOME_DIR/${SAMPLE}.sam > $SAM_TRANSCRIPTOME_DIR/stats/${SAMPLE}.stat
done

# Parse summary files because multiqc won't do it for me (how dare they)
for i in $SAM_TRANSCRIPTOME_DIR/stats/*stat
do
    SAMPLE=$(basename ${i} .stat)
    echo -e "Parsing: " ${SAMPLE}.stat
    grep ^SN $SAM_TRANSCRIPTOME_DIR/stats/${SAMPLE}.stat | cut -f 2- > $SAM_TRANSCRIPTOME_DIR/stats/${SAMPLE}_SNsection.txt
    sed -i 's/\t#.*//' $SAM_TRANSCRIPTOME_DIR/stats/${SAMPLE}_SNsection.txt
    sed -i 's/:\t/,/' $SAM_TRANSCRIPTOME_DIR/stats/${SAMPLE}_SNsection.txt
done


