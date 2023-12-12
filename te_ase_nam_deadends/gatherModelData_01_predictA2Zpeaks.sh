#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch and Emily Yi
# Contact... mbb262@cornell.edu
# Date...... 2023-11-27
# Updated... 2023-11-27
#
# Description:
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
# predict_all.sh takes argument names_list. Saves genomes from names_list to directory 
# genomes/ and saves predictions to directory preds/. If genomes/ directory already 
# contains a genome, it is not re-fetched or overwritten.

# call_peaks.sh:
# Removes scaffolds, trims 150 bp from each end of a2z predictions, filters with 0
# .75 accessibility threshold, and then does peak calling.
# takes argument names_list. Expects directory genomes/ to have the unzipped .fa 
# genomes in names_list and directory preds/ to hold a2z predictions. Requires 
# directory peaks/ to save called peaks.
# ------------------------------------------------------------------------------

# Set up conda environment -----------------------------------------------------
cd /workdir/mbb262
git clone https://github.com/twrightsman/a2z-regulatory
cd a2z-regulatory
conda env create --name a2z --file envs/a2z.yml

conda activate a2z
pip install src/python/a2z

curl 'https://zenodo.org/record/5724562/files/model-accessibility-full.h5?download=1' > model-accessibility-full.h5

# make directory for results
mkdir /workdir/mbb262/model_data/peaks

# Get data from cbsu
scp -r mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/te_ase_nam/model_data/a2z/preds /workdir/mbb262/model_data

# Save path to genomes
NAMGENOMEPATH=/workdir/mbb262/nam_hybrid_rnaseq/salmon/nam_genomes
PREDSPATH=/workdir/mbb262/model_data/preds
PEAKPATH=/workdir/mbb262/model_data/peaks
model='a2z-regulatory/model-accessibility-full.h5'


# Predict accessible regions ----------------------------------------------------

cd $NAMGENOMEPATH
for i in *.fa
do
    LONGSAMPLE=$(basename ${i} .fa)
    echo -e "Formatting: ${LONGSAMPLE}"

    # Run model
    a2z-regulatory/src/python/scripts/predict_genome.py --stride 50 $model $NAMGENOMEPATH/${LONGSAMPLE}.fa > $PREDSPATH/a2z_preds_${LONGSAMPLE}.bed
done


# Call peaks --------------------------------------------------------------------


# Loop through all genomes
cd $NAMGENOMEPATH
for i in *.fa
do
    LONGSAMPLE=$(basename ${i} .fa)
    echo -e "Formatting: ${LONGSAMPLE}"

    # Genomes are alreadu indexed so I can skip this step
    # samtools faidx ${LONGSAMPLE}.fa

    # trim 150bp off each end of the predictions and take the max prediction of all overlapping 
    # windows every 50bp and then merge touching 50bp chunks with predictions greater than or 
    # equal to 0.90 into whole regions. 
    bedtools slop -g $NAMGENOMEPATH/${LONGSAMPLE}.fa.fai -b -150 < $PREDSPATH/a2z_preds_${LONGSAMPLE}.bed > $PREDSPATH/preds_${LONGSAMPLE}.trimmed.bed
    echo -e "Done slopping off ends"

    # Call peals
    bedtools map -a <(bedtools makewindows -g $NAMGENOMEPATH/${LONGSAMPLE}.fa.fai -w 50) -b $PREDSPATH/preds_${LONGSAMPLE}.trimmed.bed -c 4 -o max | awk '$4 >= 0.75' | bedtools merge > $PEAKPATH/preds_${LONGSAMPLE}.gte90.merged.bed

done





