#!/bin/bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-07-16 
#
# Description 
#   - Script to unzip all the files in a directory using bgzip
#   - bgzip does not support globbing (i.e. bgzip *.gz)
# ---------------------------------------------------------------


for FILE in *_processed.csv.gz
do
    bgzip -d --threads 100 ${FILE}
done