# Mapping parameters ---------------
N_THREADS=40
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
ALIGN_OUT=$PROJ_DIR/output/minimap_alignments/inbred
REF_TRANS=$PROJ_DIR/references/fa_transcriptomes
FASTQ_TRIM_MERGE_DIR=$PROJ_DIR/reads/trimmedMerged
## Run minimap2 ----------------------
minimap2 -ax sr -t $N_THREADS $REF_TRANS/PHB47_v1.asm_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R002_CKDL210018333-2a-AK11819-AK30533_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R002_CKDL210018333-2a-AK11819-AK30533_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/5_Buckler-PHZ51_mecatErrorCorrected.contigs_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R003_CKDL210018333-2a-AK30535-AK30534_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R003_CKDL210018333-2a-AK30535-AK30534_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Il14H-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R004_CKDL210018333-2a-AK7741-AK25811_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R004_CKDL210018333-2a-AK7741-AK25811_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML103-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R005_CKDL210018333-2a-AK30537-AK30536_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R005_CKDL210018333-2a-AK30537-AK30536_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-NC358-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R006_CKDL210018333-2a-AK27680-AK30538_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R006_CKDL210018333-2a-AK27680-AK30538_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Ms71-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R007_CKDL210018333-2a-AK14231-AK30539_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R007_CKDL210018333-2a-AK14231-AK30539_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML247-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R008_CKDL210018333-2a-AK30540-AK16431_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R008_CKDL210018333-2a-AK30540-AK16431_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Oh43-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R009_CKDL210018333-2a-AK30542-AK30541_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R009_CKDL210018333-2a-AK30542-AK30541_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Mo17-REFERENCE-CAU-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R010_CKDL210018333-2a-7UDI301-AK1847_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R010_CKDL210018333-2a-7UDI301-AK1847_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/ Zm-HP301-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R011_CKDL210018333-2a-AK19504-5UDI1485_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R011_CKDL210018333-2a-AK19504-5UDI1485_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-B97-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R012_CKDL210018333-2a-7UDI304-AK690_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R012_CKDL210018333-2a-7UDI304-AK690_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-P39-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R013_CKDL210018333-2a-7UDI254-AK2329_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R013_CKDL210018333-2a-7UDI254-AK2329_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Mo18W-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R014_CKDL210018333-2a-AK19518-AK13670_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R014_CKDL210018333-2a-AK19518-AK13670_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Ki3-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R015_CKDL210018333-2a-AK2044-AK1952_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R015_CKDL210018333-2a-AK2044-AK1952_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Tx303-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R016_CKDL210018333-2a-GH12-GF11_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R016_CKDL210018333-2a-GH12-GF11_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-NC350-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R017_CKDL210018333-2a-AK1868-AK1852_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R017_CKDL210018333-2a-AK1868-AK1852_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-M37W-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R018_CKDL210018333-2a-AK844-AK1878_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R018_CKDL210018333-2a-AK844-AK1878_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML277-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R019_CKDL210018333-2a-GG03-7UDI222_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R019_CKDL210018333-2a-GG03-7UDI222_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Ky21-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R021_CKDL210018333-2a-AK2147-7UDI223_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R021_CKDL210018333-2a-AK2147-7UDI223_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML333-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R022_CKDL210018333-2a-AK900-7UDI293_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R022_CKDL210018333-2a-AK900-7UDI293_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML322-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R023_CKDL210018333-2a-AK2161-AK1846_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R023_CKDL210018333-2a-AK2161-AK1846_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Oh7B-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R024_CKDL210018333-2a-AK907-AK1864_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R024_CKDL210018333-2a-AK907-AK1864_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Tzi8-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R025_CKDL210018333-2a-AK903-GH07_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R025_CKDL210018333-2a-AK903-GH07_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Ki11-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R026_CKDL210018333-2a-7UDI763-AK19481_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R026_CKDL210018333-2a-7UDI763-AK19481_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML228-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R027_CKDL210018333-2a-AK18848-AK7115_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R027_CKDL210018333-2a-AK18848-AK7115_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML69-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R028_CKDL210018333-2a-AK19520-SCI5E038_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R028_CKDL210018333-2a-AK19520-SCI5E038_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-M162W-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R029_CKDL210018333-2a-AK869-AK847_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R029_CKDL210018333-2a-AK869-AK847_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Il14H-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R033_CKDL210018333-2a-AK30546-AK25858_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R033_CKDL210018333-2a-AK30546-AK25858_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Oh43-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R034_CKDL210018333-2a-AK30548-AK30547_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R034_CKDL210018333-2a-AK30548-AK30547_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Ms71-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R035_CKDL210018333-2a-AK10486-AK9660_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R035_CKDL210018333-2a-AK10486-AK9660_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Mo17-REFERENCE-CAU-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R036_CKDL210018333-2a-AK30550-AK30549_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R036_CKDL210018333-2a-AK30550-AK30549_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Ki3-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R037_CKDL210018333-2a-AK30552-AK30551_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R037_CKDL210018333-2a-AK30552-AK30551_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Tx303-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R038_CKDL210018333-2a-AK30553-AK15571_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R038_CKDL210018333-2a-AK30553-AK15571_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML247-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R039_CKDL210018333-2a-AK30554-AK27741_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R039_CKDL210018333-2a-AK30554-AK27741_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-NC350-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R040_CKDL210018333-2a-AK10434-AK30555_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R040_CKDL210018333-2a-AK10434-AK30555_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R041_CKDL210018333-2a-7UDI2717-AK23823_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R041_CKDL210018333-2a-7UDI2717-AK23823_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-M37W-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R042_CKDL210018333-2a-AK7701-AK30556_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R042_CKDL210018333-2a-AK7701-AK30556_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML277-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R043_CKDL210018333-2a-AK30558-AK30557_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R043_CKDL210018333-2a-AK30558-AK30557_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Ky21-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R044_CKDL210018333-2a-X176-AK30559_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R044_CKDL210018333-2a-X176-AK30559_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML69-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R045_CKDL210018333-2a-AK10235-AK592_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R045_CKDL210018333-2a-AK10235-AK592_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R046_CKDL210018333-2a-AK30561-AK30560_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R046_CKDL210018333-2a-AK30561-AK30560_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/PHB47_v1.asm_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R047_CKDL210018333-2a-AK30562-7UDI629_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R047_CKDL210018333-2a-AK30562-7UDI629_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/5_Buckler-PHZ51_mecatErrorCorrected.contigs_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R048_CKDL210018333-2a-K30030-AK30563_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R048_CKDL210018333-2a-K30030-AK30563_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML103-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R049_CKDL210018333-2a-AK19627-AK30564_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R049_CKDL210018333-2a-AK19627-AK30564_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/ Zm-HP301-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R051_CKDL210018333-2a-AK30567-AK30566_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R051_CKDL210018333-2a-AK30567-AK30566_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-B97-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R052_CKDL210018333-2a-AK5231-AK4415_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R052_CKDL210018333-2a-AK5231-AK4415_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-P39-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R053_CKDL210018333-2a-AK10568-AK12757_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R053_CKDL210018333-2a-AK10568-AK12757_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Mo18W-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R054_CKDL210018333-2a-AK30569-AK30568_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R054_CKDL210018333-2a-AK30569-AK30568_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R055_CKDL210018333-2a-AK30571-AK30570_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R055_CKDL210018333-2a-AK30571-AK30570_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML333-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R056_CKDL210018333-2a-AK2209-AK15971_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R056_CKDL210018333-2a-AK2209-AK15971_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML322-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R057_CKDL210018333-2a-AK30572-AK2400_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R057_CKDL210018333-2a-AK30572-AK2400_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Oh7B-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R058_CKDL210018333-2a-AK2486-AK30573_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R058_CKDL210018333-2a-AK2486-AK30573_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Tzi8-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R059_CKDL210018333-2a-AK30575-AK30574_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R059_CKDL210018333-2a-AK30575-AK30574_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Ki11-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R060_CKDL210018333-2a-AK30576-AK7506_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R060_CKDL210018333-2a-AK30576-AK7506_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML228-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R061_CKDL210018333-2a-AK30578-AK30577_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R061_CKDL210018333-2a-AK30578-AK30577_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML247-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R062_CKDL210018333-2a-AK847-AK869_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R062_CKDL210018333-2a-AK847-AK869_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-M162W-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R063_CKDL210018333-2a-AK1847-7UDI301_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R063_CKDL210018333-2a-AK1847-7UDI301_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML247-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R064_CKDL210018333-2a-AK30560-AK30561_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R064_CKDL210018333-2a-AK30560-AK30561_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/PHB47_v1.asm_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R065_CKDL210018333-2a-7UDI629-AK30562_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R065_CKDL210018333-2a-7UDI629-AK30562_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/5_Buckler-PHZ51_mecatErrorCorrected.contigs_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R066_CKDL210018333-2a-AK30563-K30030_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R066_CKDL210018333-2a-AK30563-K30030_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Oh43-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R067_CKDL210018333-2a-AK30564-AK19627_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R067_CKDL210018333-2a-AK30564-AK19627_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML103-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R068_CKDL210018333-2a-AK30565-AK20051_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R068_CKDL210018333-2a-AK30565-AK20051_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-NC358-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R069_CKDL210018333-2a-AK30566-AK30567_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R069_CKDL210018333-2a-AK30566-AK30567_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R070_CKDL210018333-2a-AK4415-AK5231_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R070_CKDL210018333-2a-AK4415-AK5231_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/ Zm-HP301-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R071_CKDL210018333-2a-AK12757-AK10568_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R071_CKDL210018333-2a-AK12757-AK10568_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-P39-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R072_CKDL210018333-2a-AK30568-AK30569_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R072_CKDL210018333-2a-AK30568-AK30569_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML52-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R073_CKDL210018333-2a-AK30570-AK30571_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R073_CKDL210018333-2a-AK30570-AK30571_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Mo18W-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R074_CKDL210018333-2a-AK15971-AK2209_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R074_CKDL210018333-2a-AK15971-AK2209_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-M37W-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R075_CKDL210018333-2a-AK2400-AK30572_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R075_CKDL210018333-2a-AK2400-AK30572_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R076_CKDL210018333-2a-AK30573-AK2486_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R076_CKDL210018333-2a-AK30573-AK2486_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML333-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R078_CKDL210018333-2a-AK7506-AK30576_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R078_CKDL210018333-2a-AK7506-AK30576_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Tzi8-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R079_CKDL210018333-2a-AK1864-AK907_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R079_CKDL210018333-2a-AK1864-AK907_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Ki11-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R080_CKDL210018333-2a-7UDI292-AK1334_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R080_CKDL210018333-2a-7UDI292-AK1334_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML69-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R081_CKDL210018333-2a-AK1867-7UDI237_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R081_CKDL210018333-2a-AK1867-7UDI237_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML247-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R082_CKDL210018333-2a-GD06-AK845_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R082_CKDL210018333-2a-GD06-AK845_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/PHB47_v1.asm_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R083_CKDL210018333-2a-AK30577-AK30578_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R083_CKDL210018333-2a-AK30577-AK30578_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/5_Buckler-PHZ51_mecatErrorCorrected.contigs_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R084_CKDL210018333-2a-AK30524-AK30523_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R084_CKDL210018333-2a-AK30524-AK30523_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Il14H-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R085_CKDL210018333-2a-AK30525-AK21985_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R085_CKDL210018333-2a-AK30525-AK21985_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Ms71-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R086_CKDL210018333-2a-AK12496-AK30526_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R086_CKDL210018333-2a-AK12496-AK30526_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Mo17-REFERENCE-CAU-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R088_CKDL210018333-2a-AK1952-AK2044_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R088_CKDL210018333-2a-AK1952-AK2044_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-B97-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R089_CKDL210018333-2a-AK985-AK987_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R089_CKDL210018333-2a-AK985-AK987_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML52-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R090_CKDL210018333-2a-AK1872-AK898_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R090_CKDL210018333-2a-AK1872-AK898_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Ki3-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R091_CKDL210018333-2a-AK1878-AK844_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R091_CKDL210018333-2a-AK1878-AK844_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Tx303-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R092_CKDL210018333-2a-AK30529-AK30528_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R092_CKDL210018333-2a-AK30529-AK30528_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-NC350-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R093_CKDL210018333-2a-AK30531-AK30530_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R093_CKDL210018333-2a-AK30531-AK30530_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML277-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R094_CKDL210018333-2a-5UDI388-AK7296_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R094_CKDL210018333-2a-5UDI388-AK7296_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML322-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R095_CKDL210018333-2a-AK16359-AK30532_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R095_CKDL210018333-2a-AK16359-AK30532_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R096_CKDL210018333-2a-GF11-GH12_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R096_CKDL210018333-2a-GF11-GH12_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-Oh7B-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R097_CKDL210018333-2a-AK1874-AK892_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R097_CKDL210018333-2a-AK1874-AK892_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-CML228-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R098_CKDL210018333-2a-AK692-GF10_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R098_CKDL210018333-2a-AK692-GF10_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 
minimap2 -ax sr -t $N_THREADS $REF_TRANS/Zm-M162W-REFERENCE-NAM-1.0_shortNames.fasta $FASTQ_TRIM_MERGE_DIR/MS21R099_CKDL210018333-2a-7UDI274-AK234_HH5V7DSX2_L1_trimMerged.fq.gz > $ALIGN_OUT/MS21R099_CKDL210018333-2a-7UDI274-AK234_HH5V7DSX2_L1_trimMerged.fq.gz_minimap2.sam 