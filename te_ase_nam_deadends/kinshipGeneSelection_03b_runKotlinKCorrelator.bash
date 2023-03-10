# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-09-26 
# Updated... 2022-09-26
#
# Description:
# Run a kotlin script from the command line
# ------------------------------------------------------------------------------

# Set up kotlin
wget https://github.com/JetBrains/kotlin/releases/download/v1.7.10/kotlin-compiler-1.7.10.zip
unzip kotlin-compiler-1.7.10.zip 
export PATH=/workdir/mbb262/kotlinc_1.7.10/bin:$PATH
export _JAVA_OPTIONS=-Xmx50g


# Run kotlin script - Terry's kinship correlator
# Generate these files within: 
kotlinc -script \
    /home/mbb262/git_projects/te_ase_nam/src/calculate_correlate_kinship.main.kts \
    /workdir/mbb262/genotypes/hmp321_282_agpv5_inbred_hybrid_all_chroms.vcf.gz \
    /workdir/mbb262/k_files/all_genes.txt \
    /workdir/mbb262/k_files/gene_coords_500bp.txt \
    /workdir/mbb262/k_files/te_gene_overlaps.txt \
    /workdir/mbb262/k_files/kinship_correlation_500bpgene_5kbte.txt

