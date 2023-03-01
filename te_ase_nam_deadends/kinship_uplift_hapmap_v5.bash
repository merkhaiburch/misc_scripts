# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-07-27 
# Updated... 2022-07-27

# Description 
# Uplift Hapmap 3.2.1 v4 SNPs to v5 using crossover package
# SNPs originate from Guillaume's heterosis paper
# ---------------------------------------------------------------

# Set up directory structure
mkdir -p /workdir/mbb262/genotypes
cd /workdir/mbb262/genotypes

# Gather INBRED Hapmap 3.2.1 data from blfs1
scp -r mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/gr226/Ramstein_AmesNAMHybrids_2019/Hmp321 /workdir/mbb262/genotypes

# Gather HYBRID Hapmap 3.2.1 data from blfs1
scp -r mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/gr226/Ramstein_AmesNAMHybrids_2019/Ames/genotypes/hybrid_imputed /workdir/mbb262/genotypes

# Gather v5 reference genome
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
gunzip Zm-B73-REFERENCE-NAM-5.0.fa.gz

# Get v4 to v5 chain file
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/chain_files/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain	


# Uplift coordinates ---------------------------------------------------------------------------------------

# Set path
export PYTHONPATH=/programs/CrossMap-0.3.8/lib64/python3.6/site-packages:/programs/CrossMap-0.3.8/lib/python3.6/site-packages
export PATH=/programs/CrossMap-0.3.8/bin:$PATH

# Create a loop that iterates through each chromosome & uplifts for INBRED data
mkdir /workdir/mbb262/genotypes/inbred_v5
for CHROM in {1..10}
do
  echo "I am on chr ${CHROM}"

  CrossMap.py vcf \
  ./B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain \
  /workdir/mbb262/genotypes/Hmp321/hmp321_282_agpv4_merged_chr${CHROM}.imputed.vcf.gz \
  ./Zm-B73-REFERENCE-NAM-5.0.fa \
  /workdir/mbb262/genotypes/inbred_v5/hmp321_282_agpv5_uplifted_merged_chr${CHROM}_imputed.vcf

  bgzip --threads 20 /workdir/mbb262/genotypes/inbred_v5/hmp321_282_agpv5_uplifted_merged_chr${CHROM}_imputed.vcf

  echo "I just finished chr ${CHROM}"
done


# Create a loop that iterates through each chromosome & uplifts for HYBRID data
mkdir /workdir/mbb262/genotypes/hybrid_v5
for CHROM in {1..10}
do
  echo "I am on chr ${CHROM}"

  CrossMap.py vcf \
  ./B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain \
  /workdir/mbb262/genotypes/hybrid_imputed/AGPv4_Ames_chr${CHROM}.PHZ51_B47-het.vcf.gz \
  ./Zm-B73-REFERENCE-NAM-5.0.fa \
  /workdir/mbb262/genotypes/hybrid_v5/AGPv5_uplifted_Ames_chr${CHROM}.PHZ51_B47-het.vcf

  bgzip --threads 20 /workdir/mbb262/genotypes/hybrid_v5/AGPv5_uplifted_Ames_chr${CHROM}.PHZ51_B47-het.vcf

  echo "I just finished chr ${CHROM}"
done


# Backup data -----------------------------------------------------------------------------------------------

# Transfer data back to blfs1
scp /workdir/mbb262/genotypes/inbred_v5/hmp321_282_agpv5_uplifted_merged_chr* mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/genotypes/Hapmap321_Ames_v5_uplifted/Ames_hybrids_v5_uplifted
scp /workdir/mbb262/genotypes/hybrid_v5/AGPv5_uplifted_Ames_chr* mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/genotypes/Hapmap321_Ames_v5_uplifted/Hapmap321_inbreds_v5_uplifted

