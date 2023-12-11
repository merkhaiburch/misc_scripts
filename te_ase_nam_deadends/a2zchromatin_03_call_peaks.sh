#!/bin/bash

# Predictions from a2z model, predictions in 600bp windows sliding 50bp at a time.
# https://github.com/twrightsman/a2z-regulatory

# ./preds contains a2z predictions
# ./peaks contains peaks called with bedtools, details in scripts

# All NAM genomes + Mo17 + PHZ51 + PHB47
# 29 genomes total

# Run setup.sh once
# If not the first time, run: conda activate a2z
#    before running prediction scripts
# Requires directory a2z-regulatory

# predict_genome.sh takes one .gz genome as argument
# predict_all.sh takes argument names_list. Saves genomes from names_list to directory genomes/ and saves predictions to directory preds/. If genomes/ directory already contains a genome, it is not re-fetched or overwritten.

# call_peaks.sh:
# Removes scaffolds, trims 150 bp from each end of a2z predictions, filters with 0.75 accessibility threshold, and then does peak calling.
# takes argument names_list. Expects directory genomes/ to have the unzipped .fa genomes in names_list and directory preds/ to hold a2z predictions. Requires directory peaks/ to save called peaks.


names_list=$1
if [ ! -f "$names_list" ]; then
    echo "$names_list does not exist."
    exit 1
fi

cat $names_list | while read name; do
    
    if [ $name == "5_Buckler-PHZ51_mecatErrorCorrected.contigs" ] || [ $name == "PHB47_v1.asm" ] || [ $name == "Zm-Mo17-REFERENCE-CAU-1.0" ]; then
	echo "Skipping $name"
	continue
    fi
    
    echo "Beginning $name"
    unzipped_genome="genomes/${name}.fa"
    if [[ ! -f $unzipped_genome ]]  # no .fa
    then
	echo "$path_to_genome does not exist"
	exit 1
    fi

    samtools faidx $unzipped_genome
    a2z_preds=preds/a2z_preds_$name.bed

    # Remove scaffolds
    sed '/^scaf/ d' < $a2z_preds > no_scaf_preds_$name.bed
    no_scaf=no_scaf_preds_$name.bed

    # trim 150 bp from each end
    bedtools slop -g $unzipped_genome.fai -b -150 < $no_scaf > preds_$name.trimmed.bed

    bedtools map -a <(bedtools makewindows -g $unzipped_genome.fai -w 50) -b preds_$name.trimmed.bed -c 4 -o max | awk '$4 >= 0.75' | bedtools merge > peaks/a2z_preds_$name.gte75.merged.bed
    
    
    rm $no_scaf
    rm preds_$name.trimmed.bed
    
    echo "Finished $name"
    done

