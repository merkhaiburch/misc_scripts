#!/bin/bash

# script assumes that kotlinc is in your PATH and that 
# the kotlin script is in the same directory as this shell script

SAM_1="simulated_reads/hybridB73Mo17_to_b73_minimap2.sam"
SAM_2="simulated_reads/hybridB73Mo17_to_mo17_minimap2.sam"

# suffixes are the reference name added to the end of the transcript name
# for each transcript/contig
SUFFIX_1="_Zm-B73-REFERENCE-NAM-5.0"
SUFFIX_2="_Zm-Mo17-REFERENCE-CAU-1.0"

# file containing a list of transcript id's to count
# e.g. Zm00001eb233020_T002
# one per line of the file
TRANSCRIPT_IDS="simulated_reads/transcript_ids.txt"

# output file
OUTPUT="out.txt"

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

# run the script
/home/mbb262/bioinformatics/kotlinc/bin/kotlinc -script count_rnaseq_reads_hybrid.main.kts -- \
	-i $SAM_1 -j $SAM_2 \
	-p $SUFFIX_1 -q $SUFFIX_2 \
	-t $TRANSCRIPT_IDS \
	-o $OUTPUT \
	-a $MIN_ALIGN \
	-l $LAX_NM -m $STRICT_NM \
	-n $MIN_MAPQ