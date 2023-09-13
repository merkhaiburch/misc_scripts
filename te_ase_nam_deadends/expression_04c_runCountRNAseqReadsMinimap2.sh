# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-10-03
# Updated... 2023-02-28
#
# Description:
# Run count_rnaseq_reads_minimap2.kt
# ------------------------------------------------------------------------------

# Set up kotlin
# wget https://github.com/JetBrains/kotlin/releases/download/v1.7.10/kotlin-compiler-1.7.10.zip
# unzip kotlin-compiler-1.7.10.zip 
# export PATH=/workdir/mbb262/kotlinc_1.7.10/bin:$PATH
# export _JAVA_OPTIONS=-Xmx50g


# Filter inbred alignments to primary alignments
SAM_TRANSCRIPTOME_DIR=/workdir/mbb262/nam_hybrid_rnaseq/output/minimap_alignments/inbred

mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/output/minimap_alignments/inbred_primary
PRIMARY_SAM_DIR=/workdir/mbb262/nam_hybrid_rnaseq/output/minimap_alignments/inbred_primary

for i in $SAM_TRANSCRIPTOME_DIR/*.sam
do
    SAMPLE=$(basename ${i} .sam)

    # Get primary alignments from sam
    /programs/samtools-1.15.1-r/bin/samtools view -e 'flag <= 16' $SAM_TRANSCRIPTOME_DIR/${SAMPLE}.sam -o $PRIMARY_SAM_DIR/${SAMPLE}_primary.sam
done

# Run Kotlin script on inbreds - from Zack
/home/mbb262/bioinformatics/kotlinc/bin/kotlinc -script \
    /home/mbb262/git_projects/te_ase_nam/src/04b_count_rnaseq_reads_minimap2.main.kts


# Run Kotlin script on hybreds - from Ana
/home/mbb262/bioinformatics/kotlinc/bin/kotlinc -script \
    /home/mbb262/git_projects/te_ase_nam/src/TODO

