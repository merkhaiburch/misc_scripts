#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-12-11
# Updated... 2022-12-11
#
# Description:
# Subsample sites in a genotype file, export the site summary from TASSEL
# ------------------------------------------------------------------------------


for FILE in *.gz
do
	/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/output_counts/goodman_filtered_count/debug.txt \
    -Xmx200g \
    -maxThreads 60 \
    -importGuess ./$FILE \
    -noDepth \
    -filterAlign \
    -subsetSites 400000 \
    -genotypeSummary site \
    -export ${FILE}_subsample_sites.txt
done