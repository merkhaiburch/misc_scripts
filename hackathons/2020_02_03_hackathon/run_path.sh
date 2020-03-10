#!/bin/bash

# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-02-18 
#
# Description 
#   - Script by Evan Long to 
#   - 
# ---------------------------------------------------------------

# Notes
# $1 config
# $2 consensus method
# $3 fastq Dir

# Provides NAM haplotypes after running consensus from DB, methods file maxDiv = 0.001 
/workdir/eml255/tools/tassel-5-standalone/run_pipeline.pl \
  -Xmx100g \
  -debug \
  -HaplotypeGraphBuilderPlugin \
  -configFile $1 \
  -methods $2 \
  -includeVariantContexts false \
  -endPlugin \
  -WriteFastaFromGraphPlugin \
  -outputFile phg_pangenome \
  -endPlugin > PHG_pangeome.log


# Indexes haplotype fastqs using minimap
/programs/minimap2-2.17/minimap2 \
  -d pangenome_index  \
  -k 21 \
  -w 11 \
  -I 90G phg_pangenome.fa


# Mapping 282 fastqs onto pangenome haplotypes (mapping 282 reads to collapsed PHG haplotypes)
/workdir/eml255/tools/tassel-5-standalone/run_pipeline.pl \
  -debug \
  -Xmx200g \
  -configParameters $1 \
  -HaplotypeGraphBuilderPlugin \
  -configFile $1 \
  -methods $2 \
  -includeVariantContexts false \
  -includeSequences false  \
  -endPlugin \
  -FastqToMappingPlugin \
  -minimap2IndexFile pangenome_index \
  -minimapLocation /programs/minimap2-2.17/minimap2 \
  -keyFile keyfile.txt \
  -fastqDir $3 \
  -methodName "dummy2" \
  -methodDescription "dummy2" \
  -debugDir debug/ \
  -endplugin > PHG_dipPath.log


#Find Paths
/workdir/eml255/tools/tassel-5-standalone/run_pipeline.pl \
  -debug \
  -Xmx240g \
  -configParameters $1 \
  -HaplotypeGraphBuilderPlugin \
  -configFile $1 \
  -methods $2,refRegionGroup \
  -includeVariantContexts false \
  -includeSequences false \
  -endPlugin \
  -BestHaplotypePathPlugin \
  -keyFile keyfile_pathing.txt \
  -outDir out_path \
  -readMethod "dummy2" \
  -pathMethod "dummy2" \
  -endPlugin > PHG_Path.log


# export paths
/workdir/tassel-5-standalone/run_pipeline.pl  \
  -debug \
  -Xmx100g \
  -configParameters config.txt \
  -HaplotypeGraphBuilderPlugin \
  -configFile config.txt \
  -methods mummer4,refRegionGroup \
  -includeVariantContexts false \
  -includeSequences false -endPlugin \
  -BestHaplotypePathPlugin \
  -keyFile keyfile_pathing.txt \
  -outDir out_path \
  -readMethod "dummy2" \
  -pathMethod "dummy2" \
  -endPlugin








