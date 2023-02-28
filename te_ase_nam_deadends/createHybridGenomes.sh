#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-02-23
# Updated... 2023-02-23
#
# Description:
# Create "artificial" hybrids by merging genomes together and adding source genome ID
# ------------------------------------------------------------------------------

# Create hybrid genomes ------------------------------------------

# Hybrids are: NAM x B73, NAM x PHZ51, NAM x PHB47, PHB47 x PHZ51
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/references/combined_fa_transcriptomes
B73=$FA_TRANSCRIPTOME_DIR/Zm-B73-REFERENCE-NAM-5.0_B73v5_mRNA_axSplice_eqx_canonical_named.fa
PHZ51=$FA_TRANSCRIPTOME_DIR/5_Buckler-PHZ51_mecatErrorCorrected.contigs.fasta_B73v5_mRNA_axSplice_eqx_canonical_named.fa
PHB47=$FA_TRANSCRIPTOME_DIR/PHB47_gapFilled.fasta_B73v5_mRNA_axSplice_eqx_canonical_named.fa
FA_TRANSCRIPTOME_DIR=/workdir/mbb262/nam_hybrid_rnaseq/references/fa_transcriptomes
FA_TRANSCRIPTOME_COMBINED_DIR=/workdir/mbb262/nam_hybrid_rnaseq/references/combined_fa_transcriptomes

# SAMPLE=$(basename ${i} .fa)

# Create all possible hybrids with these three founders
for i in $FA_TRANSCRIPTOME_DIR/*.fa
do    
    SAMPLE=$(basename ${i} .fa)
    echo "Create super genome for: ${SAMPLE}"

    # ALL x B73
    cat ${B73} $FA_TRANSCRIPTOME_DIR/${SAMPLE}.fa > $FA_TRANSCRIPTOME_COMBINED_DIR/${SAMPLE}.faZm-B73-REFERENCE-NAM-5.0_B73v5_mRNA_axSplice_eqx_canonical_named.fa

    # ALL x PHZ51
    cat ${PHZ51} $FA_TRANSCRIPTOME_DIR/${SAMPLE}.fa > $FA_TRANSCRIPTOME_COMBINED_DIR/${SAMPLE}.fa5_Buckler-PHZ51_mecatErrorCorrected.contigs.fasta_B73v5_mRNA_axSplice_eqx_canonical_named.fa

    # ALL x PHB47
    cat ${PHB47} \
        $FA_TRANSCRIPTOME_DIR/${SAMPLE}.fa > $FA_TRANSCRIPTOME_COMBINED_DIR/${SAMPLE}.faPHB47_gapFilled.fasta_B73v5_mRNA_axSplice_eqx_canonical_named.fa
done

# Copy inbred transcriptomes to combined transcriptomes
cp $FA_TRANSCRIPTOME_DIR/* $FA_TRANSCRIPTOME_COMBINED_DIR/

# Remove duplicates
cd /workdir/mbb262/nam_hybrid_rnaseq/references/combined_fa_transcriptomes
rm 5_Buckler-PHZ51_mecatErrorCorrected.contigs.fasta_B73v5_mRNA_axSplice_eqx_canonical_named.fa5_Buckler-PHZ51_mecatErrorCorrected.contigs.fasta_B73v5_mRNA_axSplice_eqx_canonical_named.fa
rm PHB47_gapFilled.fasta_B73v5_mRNA_axSplice_eqx_canonical_named.faPHB47_gapFilled.fasta_B73v5_mRNA_axSplice_eqx_canonical_named.fa
rm Zm-B73-REFERENCE-NAM-5.0_B73v5_mRNA_axSplice_eqx_canonical_named.faZm-B73-REFERENCE-NAM-5.0_B73v5_mRNA_axSplice_eqx_canonical_named.fa
ll | wc -l # expecting 113 files
