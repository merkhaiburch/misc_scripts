#!/bin/bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-06-14 
#
# Description 
#   - kmer distribution group for Hackathon
#   - trim 10 PanAnd species reads using Trimomatic and check 
#	- quality using FastQC
# ---------------------------------------------------------------


# -------------------------------------
# Use trimmotatic
# -------------------------------------

# Zea mays
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1=zea_mays_1.fq.gz \
    in2=zea_mays_2.fq.gz \
    out1=zea_mays_clean1.fq \
    out2=zea_mays_clean2.fq

# Zea diploperennis (Leaf from Kellogg PI 441930) [?]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1=zea_diploperennis_1.fq.gz \
    in2=zea_diploperennis_2.fq.gz \
    out1=zea_diploperennis_1_clean1.fq \
    out2=zea_diploperennis_2_clean2.fq

# Tripsacum dactyloides var floridanum (Buckler lab clone T007-0002) [?]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1=tripsacum_dactyloides_floridanum_1.fq.gz \
    in2=tripsacum_dactyloides_floridanum_2.fq.gz \
    out1=tripsacum_dactyloides_floridanum_1_clean1.fq \
    out2=tripsacum_dactyloides_floridanum_2_clean2.fq

# Vossia cuspidata (Buckler lab clone A1077-0002) [?]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1=vossia_cuspidata_1.fq.gz \
    in2=vossia_cuspidata_2.fq.gz \
    out1=vossia_cuspidata_1_clean1.fq \
    out2=vossia_cuspidata_2_clean2.fq

# Chrysopogon aciculatus (AUB 135) [AN21TSTL0006]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1=chrysopogon_aciculatus_1.fq.gz \
    in2=chrysopogon_aciculatus_2.fq.gz \
    out1=chrysopogon_aciculatus_1_clean1.fq \
    out2=chrysopogon_aciculatus_2_clean2.fq

# Themeda avenacea (AQ0476980) [AN20TSCR000216] --> next one to download
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1=themeda_avenacea_1.fq.gz \
    in2=themeda_avenacea_2.fq.gz \
    out1=themeda_avenacea_1_clean1.fq \
    out2=themeda_avenacea_2_clean2.fq

# Andropogon gerardii (USA: Kansas, Ellis) [AN20NSCR000358]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1=andropogon_gerardii_1.fq.gz \
    in2=andropogon_gerardii_2.fq.gz \
    out1=andropogon_gerardii_1_clean1.fq \
    out2=andropogon_gerardii_2_clean2.fq

# Ischaemum koenigii (Pasquet 1196) [AN21TNTL0138]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1=ischaemum_koenigii_1.fq.gz \
    in2=ischaemum_koenigii_2.fq.gz \
    out1=ischaemum_koenigii_1_clean1.fq \
    out2=ischaemum_koenigii_2_clean2.fq

# Coix lacryma-jobi (PI 324509) [AN21TSTL0025]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1=coix_lacryma_jobi_1.fq.gz \
    in2=coix_lacryma_jobi_2.fq.gz \
    out1=coix_lacryma_jobi_1_clean1.fq \
    out2=coix_lacryma_jobi_2_clean2.fq

# Arthraxon junnarensis (Kew 88886) [AN21TSTL0030]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1=arthraxon_junnarensis_1.fq.gz \
    in2=arthraxon_junnarensis_2.fq.gz \
    out1=arthraxon_junnarensis_1_clean1.fq \
    out2=arthraxon_junnarensis_2_clean2.fq

# Skipping other team's genomes


# ----------------------------------------
# Run fastqc on two random genomes
# to see if there are any enriched kmers
# ----------------------------------------

/programs/FastQC-0.11.8/fastqc ischaemum_koenigii_1.fq.gz ischaemum_koenigii_2.fq.gz

/programs/FastQC-0.11.8/fastqc zea_diploperennis_1.fq.gz zea_diploperennis_2.fq.gz

/programs/FastQC-0.11.8/fastqc schizachyrium_fragile_1.fq.gz schizachyrium_fragile_2.fq.gz
