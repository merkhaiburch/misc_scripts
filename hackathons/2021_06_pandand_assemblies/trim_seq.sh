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
    out1=_clean1.fq \
    out2=_clean2.fq

mv MCRTL020_CKDL200164854-1a-AK6510-AK6670_HMGJWDSXY_L4_1.fq.gz zea_diploperennis_1.fq.gz
mv MCRTL020_CKDL200164854-1a-AK6510-AK6670_HMGJWDSXY_L4_2.fq.gz zea_diploperennis_2.fq.gz

# Tripsacum dactyloides var floridanum (Buckler lab clone T007-0002) [?]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1= \
    in2=\
    out1=_clean1.fq \
    out2=_clean2.fq
mv T0007_USPD16080542-1_H3325CCXY_L3_1.fq.gz tripsacum_dactyloides_floridanum_1.fq.gz
mv T0007_USPD16080542-1_H3325CCXY_L3_2.fq.gz tripsacum_dactyloides_floridanum_2.fq.gz

# Vossia cuspidata (Buckler lab clone A1077-0002) [?]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1= \
    in2=\
    out1=_clean1.fq \
    out2=_clean2.fq
mv A1077002_USPD16097230-AK4973-AK12977_HJM33DSXX_L4_1.fq.gz vossia_cuspidata_1.fq.gz
mv A1077002_USPD16097230-AK4973-AK12977_HJM33DSXX_L4_2.fq.gz vossia_cuspidata_2.fq.gz

# Chrysopogon aciculatus (AUB 135) [AN21TSTL0006]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1= \
    in2=\
    out1=_clean1.fq \
    out2=_clean2.fq
mv MCRTL042_CKDL200164855-1B-AK6514-AK6680_HMJV7DSXY_L3_1.fq.gz chrysopogon_aciculatus_1.fq.gz
mv MCRTL042_CKDL200164855-1B-AK6514-AK6680_HMJV7DSXY_L3_2.fq.gz chrysopogon_aciculatus_2.fq.gz

# Themeda avenacea (AQ0476980) [AN20TSCR000216] --> next one to download
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1= \
    in2=\
    out1=_clean1.fq \
    out2=_clean2.fq
mv AN20T012_CKDL200169515-1a-AK5845-AK6676_HLTF5DSXY_L3_1.fq.gz themeda_avenacea_1.fq.gz
mv AN20T012_CKDL200169515-1a-AK5845-AK6676_HLTF5DSXY_L3_2.fq.gz themeda_avenacea_2.fq.gz

# Andropogon gerardii (USA: Kansas, Ellis) [AN20NSCR000358]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1= \
    in2=\
    out1=_clean1.fq \
    out2=_clean2.fq
mv AN20N001_CKDL200149282-1a-AK4958-7UDI938_HNLFVDSXX_L1_1.fq.gz andropogon_gerardii_1.fq.gz
mv AN20N001_CKDL200149282-1a-AK4958-7UDI938_HNLFVDSXX_L1_2.fq.gz andropogon_gerardii_2.fq.gz

# Ischaemum koenigii (Pasquet 1196) [AN21TNTL0138]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1= \
    in2=\
    out1=_clean1.fq \
    out2=_clean2.fq
mv AN21TNTL0138_CKDL210009738-1a-AK6510-AK6670_H5C7LDSX2_L4_1.fq.gz ischaemum_koenigii_1.fq.gz
mv AN21TNTL0138_CKDL210009738-1a-AK6510-AK6670_H5C7LDSX2_L4_2.fq.gz ischaemum_koenigii_2.fq.gz

# Coix lacryma-jobi (PI 324509) [AN21TSTL0025]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1= \
    in2=\
    out1=_clean1.fq \
    out2=_clean2.fq

mv AN21TS25_CKDL210002323-1a-AK6532-AK6712_HVGWVDSXY_L3_1.fq.gz coix_lacryma_jobi_1.fq.gz
mv AN21TS25_CKDL210002323-1a-AK6532-AK6712_HVGWVDSXY_L3_2.fq.gz coix_lacryma_jobi_2.fq.gz

# Arthraxon junnarensis (Kew 88886) [AN21TSTL0030]
/home/mbb262/bioinformatics/bbmap/bbduk.sh \
    in1= \
    in2=\
    out1=_clean1.fq \
    out2=_clean2.fq

mv AN21TS30_CKDL210002324-1a-AK6516-AK4469_HVGWVDSXY_L4_1.fq.gz arthraxon_junnarensis_1.fq.gz
mv AN21TS30_CKDL210002324-1a-AK6516-AK4469_HVGWVDSXY_L4_2.fq.gz arthraxon_junnarensis_2.fq.gz

