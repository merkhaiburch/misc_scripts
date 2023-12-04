#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-12-01
# Updated... 2023-12-01
#
# Description:
# Code to create a blastdb and search for sequences against it
# ------------------------------------------------------------------------------

# Change directory
cd /path/to/your/working/dir

# Download genome from irods
iget /ibl/home/assemblies/andropogoneae/private/final_8_19_22/Zd-Gigi-REFERENCE-PanAnd-1.0.fasta.gz

# Unzip genome
gunzip Zd-Gigi-REFERENCE-PanAnd-1.0.fasta.gz

# Create a blast database
makeblastdb \
      -in Zd-Gigi-REFERENCE-PanAnd-1.0.fasta \
      -input_type fasta \
      -dbtype nucl \
      -title gigiRef \
      -out gigi.db

# Run blast
blastn \
      -db gigi.db \
      -query your_fasta_file_wth_all_nitrogen_genes.fa \
      -out blast_matches_tabular.txt \
      -outfmt "6 qseqid qgi qacc qdescr ssequid sacc sallseqid evalue bitscore length pident nident"


# Run minimap
/programs/minimap2-2.26/minimap2 -a \
      Zd-Gigi-REFERENCE-PanAnd-1.0.fasta \
      your_fasta_file_wth_all_nitrogen_genes.fa > alignment.sam