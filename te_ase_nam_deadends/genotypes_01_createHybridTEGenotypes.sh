# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-10-03
# Updated... 2022-10-05
#
# Description:
# Format TE genotype table to have hybrid genotypes
# ------------------------------------------------------------------------------

# Get file from blsf1
scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mcs368/anchorwave_indel/NAM_anchorwave_indels.PAV01.noAncestralSpecification.2022-09-27.txt.gz /workdir/mbb262/genotypes

# Format the above genotype table in R using: model_fit_expression_te_k.R


# Load genotypes into tassel, create hybrid genotypes --------------------------

# Variables
GENO_DIR=/workdir/mbb262/genotypes

# Run tassel to create hybrid genotypes for:
# B73 x NAM inbreds, PHZ51 x NAM inbreds, and PHB47 x NAM inbreds
/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug $GENO_DIR/debug_hybrid_TE.log \
    -Xmx200g \
    -maxThreads 30 \
    -importGuess $GENO_DIR/NAM_anchorwave_indels.PAV01.noAncestralSpecification.2022-09-27_rformat.txt \
    -CreateHybridGenotypesPlugin \
    -numericGenotype true \
    -hybridFile /workdir/mbb262/genotypes/te_ase_hybrids_for_tassel_b73_only.txt \
    -hybridChar "_" \
    -endPlugin \
    -export $GENO_DIR/hybrid_numeric_geno.txt 


# Remove the first line from these files
sed '1d' hybrid_numeric_geno.txt > hybrid_numeric_geno_no_header.txt
