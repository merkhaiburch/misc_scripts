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


git clone https://github.com/twrightsman/a2z-regulatory
cd a2z-regulatory
conda create --name a2z python
mamba env update -n a2z --file envs/a2z.yml

conda activate a2z
pip install src/python/a2z

curl 'https://zenodo.org/record/5724562/files/model-accessibility-full.h5?download=1' > model-accessibility-full.h5
