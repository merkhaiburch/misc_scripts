
# make directory for results
mkdir /workdir/mbb262/model_data/peaks

# Get data from cbsu
scp -r mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/te_ase_nam/model_data/a2z/preds /workdir/mbb262/model_data

# Save path to genomes
NAMGENOMEPATH=/workdir/mbb262/nam_hybrid_rnaseq/salmon/nam_genomes
PREDSPATH=/workdir/mbb262/model_data/preds
PEAKPATH=/workdir/mbb262/model_data/peaks

# Loop through all genomes
cd $NAMGENOMEPATH
for i in *.fa
do
    LONGSAMPLE=$(basename ${i} .fa)
    echo -e "Formatting: ${LONGSAMPLE}"

    # Genomes are alreadu indexed so I can skip this step
    # samtools faidx ${LONGSAMPLE}.fa

    # trim 150bp off each end of the predictions and take the max prediction of all overlapping 
    # windows every 50bp and then merge touching 50bp chunks with predictions greater than or 
    # equal to 0.90 into whole regions. 
    bedtools slop -g $NAMGENOMEPATH/${LONGSAMPLE}.fa.fai -b -150 < $PREDSPATH/a2z_preds_${LONGSAMPLE}.bed > $PREDSPATH/preds_${LONGSAMPLE}.trimmed.bed
    echo -e "Done slopping off ends"

    # Call peals
    bedtools map -a <(bedtools makewindows -g $NAMGENOMEPATH/${LONGSAMPLE}.fa.fai -w 50) -b $PREDSPATH/preds_${LONGSAMPLE}.trimmed.bed -c 4 -o max | awk '$4 >= 0.75' | bedtools merge > $PEAKPATH/preds_${LONGSAMPLE}.gte90.merged.bed

done