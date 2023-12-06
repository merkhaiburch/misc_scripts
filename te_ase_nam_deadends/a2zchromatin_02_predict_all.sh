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

model='a2z-regulatory/model-accessibility-full.h5'
if [ ! -f "$model" ]; then
    echo "$model does not exist."
    exit 1
fi

names_list=$1
if [ ! -f "$names_list" ]; then
    echo "$names_list does not exist."
    exit 1
fi

cat $names_list | while read name; do
    echo "Beginning $name"
    path_to_genome="genomes/${name}.fa"
    if [[ ! -f $path_to_genome ]]  # no .fa
    then
	echo $path_to_genome.gz
	if [[ ! -f $path_to_genome.gz ]]  # no .fa.gz
	then
	    wget -P genomes/ https://ftp.maizegdb.org/MaizeGDB/FTP/${name}/${name}.fa.gz
	fi
	gzip -d genomes/${name}.fa.gz
    # else, have the .fa, do nothing
    fi
    
    echo "Unzipped $name.fa.gz"
    unzipped=genomes/${name}.fa # add .fa extension
    
    a2z-regulatory/src/python/scripts/predict_genome.py --stride 50 $model $unzipped > preds/a2z_preds_$name.bed    
    echo "Predicted $name"
done
