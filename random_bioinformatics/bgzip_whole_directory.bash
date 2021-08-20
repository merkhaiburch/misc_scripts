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

# unzip
for FILE in *_processed.csv.gz
do
    bgzip -d --threads 100 ${FILE}
done


# zip
for FILE in *.txt
do
    bgzip --threads 30 ${FILE}
done

# index file (for bcftools) using bcftools
for FILE in *.vcf.gz
do
    bcftools index --threads 20 ${FILE}
done