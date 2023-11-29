# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-11-22
# Updated... 2023-11-22
#
# Description:
# Test pipeline to subset fastas to the same set of ids before mapping with salmon
# ------------------------------------------------------------------------------

# Mark which are the test files
# Subset to debug --> c("MS21R193", "MS21R254", "MS21R317", "MS21R384")
# CML322/PHB47

# Make a test directory
mkdir /workdir/mbb262/test_intersect_transcripts
mkdir /workdir/mbb262/test_intersect_transcripts/salmon
mkdir /workdir/mbb262/test_intersect_transcripts/output

TESTDIR=/workdir/mbb262/test_intersect_transcripts
SALMONDIR=/workdir/mbb262/test_intersect_transcripts/salmon
SALMONOUT=/workdir/mbb262/test_intersect_transcripts/output

# Copy over the fastas we need
cp /workdir/mbb262/nam_hybrid_rnaseq/salmon/nam_genomes/liftoff_inbred/agat_sequences/Zm-CML322-REFERENCE-NAM-1.0_agat.fa $TESTDIR
cp /workdir/mbb262/nam_hybrid_rnaseq/salmon/nam_genomes/liftoff_inbred/agat_sequences/PHB47_v1.asm_agat.fa $TESTDIR
cp /workdir/mbb262/nam_hybrid_rnaseq/salmon/nam_genomes/PHB47_v1.asm.fa $TESTDIR
cp /workdir/mbb262/nam_hybrid_rnaseq/salmon/nam_genomes/Zm-CML322-REFERENCE-NAM-1.0.fa $TESTDIR

# bash: Find fastas that agat extract sequences found, export names of all reads, save to file
cd $TESTDIR
awk 'sub(/^>/, "")' Zm-CML322-REFERENCE-NAM-1.0_agat.fa > cml322.txt
awk 'sub(/^>/, "")' PHB47_v1.asm_agat.fa > phb47.txt

# # R: Read all files into R, make into a list
# cml322 <- read.delim("/workdir/mbb262/test_intersect_transcripts/cml322.txt", header = F)
# phb47 <- read.delim("/workdir/mbb262/test_intersect_transcripts/phb47.txt", header = F)
# # R: Intersect transcript names
# lala <- data.frame(intersect(cml322$V1, phb47$V1))
# # R: Check if lengths changed
# dim(cml322)
# dim(phb47)
# length(lala)
# # R: Run another check to make sure they're present in more than 1 genome
# temp <- rbind(cml322, phb47)
# countTrans <- table(temp$V1)
# countLala <- table(lala) # should all = 1
# # R: Export common set of ids
# write.table(lala, "/workdir/mbb262/test_intersect_transcripts/cmlphbshared.txt",
#             quote = FALSE, col.names = FALSE, row.names = FALSE)

# Results in 68,714 transcripts

# bash: subset all fastas with the common set of ids
/programs/samtools-1.15.1-r/bin/samtools faidx \
    ./Zm-CML322-REFERENCE-NAM-1.0_agat.fa \
    -r ./cmlphbshared.txt \
    -o ./Zm-CML322-REFERENCE-NAM-1.0_agat_sub.fa

/programs/samtools-1.15.1-r/bin/samtools faidx \
    ./PHB47_v1.asm_agat.fa \
    -r ./cmlphbshared.txt \
    -o ./PHB47_v1.asm_agat_sub.fa


# Set up for indexing
cd $TESTDIR
for i in *1.0.fa
do
    LONGSAMPLE=$(basename ${i} .fa)
    SAMPLE=$(basename ${i} "-REFERENCE-NAM-1.0.fa" | sed 's/^Zm-//' | sed 's/^5_Buckler-//' | sed 's/_mecatErrorCorrected.contigs.fa$//' | sed 's/-REFERENCE-NAM-5.0.fa$//' | sed 's/_v1.asm.fa$//' | sed 's/_v1.asm-1.0.fa$//' | sed 's/-REFERENCE-CAU-2.0.fa$//')
    echo -e "Formatting: ${LONGSAMPLE}"
    echo -e "Short name: ${SAMPLE}"

    # Remove everything after the space
    sed 's/\s.*$//' $TESTDIR/${LONGSAMPLE}.fa > $SALMONDIR/${LONGSAMPLE}_labeled.fa

    # Add which genome it's from
    sed "/^>/s/$/&_${SAMPLE}/" -i $SALMONDIR/${LONGSAMPLE}_labeled.fa
done

# Collect chromosome and scaffold names in genomes for decoys
for i in $SALMONDIR/*_labeled.fa
do
    LONGSAMPLE=$(basename ${i} _labeled.fa)
    echo -e "Formatting: ${LONGSAMPLE}"

    grep '^>' $SALMONDIR/${LONGSAMPLE}_labeled.fa | sed 's/^>//g' > $SALMONDIR/${LONGSAMPLE}_decoys.txt
done

# Add on which genome the transcript sequences are from
for i in $TESTDIR/*_sub.fa
do
    LONGSAMPLE=$(basename ${i} _agat_sub.fa)
    SAMPLE=$(basename ${i} "-REFERENCE-NAM-1.0_agat.fa" | sed 's/^Zm-//' | sed 's/^5_Buckler-//' | sed 's/_mecatErrorCorrected.contigs_agat.fa$//' | sed 's/-REFERENCE-NAM-5.0_agat.fa$//' | sed 's/_v1.asm_agat.fa$//' | sed 's/_v1.asm_agat_sub.fa$//' | sed 's/-REFERENCE-CAU-2.0_agat.fa$//')
    echo -e "Formatting: ${LONGSAMPLE}"
    echo -e "Short name: ${SAMPLE}"

    sed "/^>/s/$/&_${SAMPLE}/" $TESTDIR/${LONGSAMPLE}_agat_sub.fa > $SALMONDIR/${LONGSAMPLE}_labeledAgatTranscripts.fa
done


# Run indexing
cat $SALMONDIR/Zm-CML322-REFERENCE-NAM-1.0_decoys.txt $SALMONDIR/PHB47_v1.asm-1.0_decoys.txt > $SALMONDIR/Zm-CML322-REFERENCE-NAM-1.0_PHB47_v1.asm-1.0.txt 

cat $SALMONDIR/Zm-CML322-REFERENCE-NAM-1.0_labeledAgatTranscripts.fa $SALMONDIR/PHB47_v1.asm_labeledAgatTranscripts.fa $SALMONDIR/Zm-CML322-REFERENCE-NAM-1.0_labeled.fa $SALMONDIR/PHB47_v1.asm-1.0_labeled.fa > $SALMONDIR/Zm-CML322-REFERENCE-NAM-1.0_PHB47_v1.asm-1.0_gentrome.fa 
conda activate salmon
salmon index --threads 40 --kmerLen 31 --transcripts $SALMONDIR/Zm-CML322-REFERENCE-NAM-1.0_PHB47_v1.asm-1.0_gentrome.fa --index $SALMONDIR/Zm-CML322-REFERENCE-NAM-1.0_PHB47_v1.asm-1.0 --decoys $SALMONDIR/Zm-CML322-REFERENCE-NAM-1.0_PHB47_v1.asm-1.0.txt --keepDuplicates 


# Run quantification
# Subset to debug --> c("MS21R193", "MS21R254", "MS21R317", "MS21R384")
FASTQ_TRIM_DIR=/workdir/mbb262/nam_hybrid_rnaseq/reads/trimmed
SALMONOUT=/workdir/mbb262/test_intersect_transcripts/output
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R193_CKDL210018333-2a-AK1880-GB07_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R193_CKDL210018333-2a-AK1880-GB07_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMONDIR/Zm-CML322-REFERENCE-NAM-1.0_PHB47_v1.asm-1.0 --validateMappings --numBootstraps 30 --output $SALMONOUT/MS21R193 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R254_CKDL210018333-2a-5UDI2095-AK19501_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R254_CKDL210018333-2a-5UDI2095-AK19501_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMONDIR/Zm-CML322-REFERENCE-NAM-1.0_PHB47_v1.asm-1.0 --validateMappings --numBootstraps 30 --output $SALMONOUT/MS21R254 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R317_CKDL210018333-2a-AK30483-AK30482_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R317_CKDL210018333-2a-AK30483-AK30482_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMONDIR/Zm-CML322-REFERENCE-NAM-1.0_PHB47_v1.asm-1.0 --validateMappings --numBootstraps 30 --output $SALMONOUT/MS21R317 
salmon quant --threads $N_THREADS --libType A -1 $FASTQ_TRIM_DIR/MS21R384_CKDL210018333-2a-AK30528-AK30529_HH5V7DSX2_L1_trimmedQC_1.fq.gz -2 $FASTQ_TRIM_DIR/MS21R384_CKDL210018333-2a-AK30528-AK30529_HH5V7DSX2_L1_trimmedQC_2.fq.gz --index $SALMONDIR/Zm-CML322-REFERENCE-NAM-1.0_PHB47_v1.asm-1.0 --validateMappings --numBootstraps 30 --output $SALMONOUT/MS21R384 

# Check stuff in R















