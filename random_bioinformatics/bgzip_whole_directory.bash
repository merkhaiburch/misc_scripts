#!/bin/bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-07-16 
# Updated... 2022-12-28
#
# Description 
#   - Script to unzip all the files in a directory using bgzip
#   - bgzip does not support globbing (i.e. bgzip *.gz)
# ---------------------------------------------------------------

# Parallel command to zip
find ./*.txt -type f > zippable_files.list
/programs/parallel/bin/parallel -j 80 "gzip {}" :::: ./zippable_files.list

# Parallel command to unzip
find ./*.gz -type f > unzippable_files.list
/programs/parallel/bin/parallel -j 11 "gunzip {}" :::: ./unzippable_files.list


# unzip
for FILE in *.txt
do
    bgzip -d --threads 20 ${FILE}
done


# zip txt file
for FILE in *.txt
do
    bgzip --threads 40 ${FILE}
done


# zip csv file
for FILE in *.csv
do
    bgzip --threads 20 ${FILE}
done


# zip vcf file
for FILE in *.vcf
do
    bgzip --threads 30 ${FILE}
done


# index file (for bcftools) using bcftools
for FILE in *.vcf.gz
do
    bcftools index --threads 20 ${FILE}
done