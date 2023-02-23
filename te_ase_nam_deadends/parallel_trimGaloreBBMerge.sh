## Trim Galore! for paired end reads -------------------------------------------------------------

# Make directories
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/reads/raw
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/reads/trimmed
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/reads/trimming_reports
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/output/fastqc
mkdir -p /workdir/mbb262/nam_hybrid_rnaseq/reads/parallel_trimmed

# To prove that parallel works with these paired reads
# https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/53688-gnu-parallel-cutadapt-with-paired-end-reads
find *_L1_1.fq.gz | sed 's/_1.fq.gz$//' | parallel 'echo {}_1.fq.gz {}_2.fq.gz'

# Parallel version, much faster, DRY RUN
# /programs/parallel/bin/parallel --dry-run -j 10 "trim_galore --cores 3 --paired ./{}_1.fq.gz ./{}_2.fq.gz --length 148 --output_dir ./ --quality 20 --fastqc --fastqc_args '-o ./' --basename {}" :::: <(find *_L1_1.fq.gz | sed 's/_1.fq.gz$//')

# Actual parallel version
/programs/parallel/bin/parallel -j 12 "trim_galore --cores 5 --paired ./{}_1.fq.gz ./{}_2.fq.gz --length 60 --output_dir /workdir/mbb262/nam_hybrid_rnaseq/reads/trimmed/ --fastqc --fastqc_args '-o ./' --basename {}" :::: <(find *_L1_1.fq.gz | sed 's/_1.fq.gz$//')


# move summary outputs to a differnt directory
mv /workdir/mbb262/nam_hybrid_rnaseq/reads/trimmed/*_trimming_report.txt /workdir/mbb262/nam_hybrid_rnaseq/reads/trimming_reports
mv /workdir/mbb262/nam_hybrid_rnaseq/reads/raw/*_fastqc* /workdir/mbb262/nam_hybrid_rnaseq/reads/trimming_reports



##  Merge paired reads with BBmerge -------------------------------------------------------------------

# Variables
mkdir /workdir/mbb262/nam_hybrid_rnaseq/reads/merged_reports
PROJ_DIR=/workdir/mbb262/nam_hybrid_rnaseq
FASTQ_TRIM=$PROJ_DIR/reads/trimmed
MERGE_DIR=$PROJ_DIR/reads/merged
OUT_DIR_INSERT_SIZE=$PROJ_DIR/reads/merged_reports

for i in $FASTQ_TRIM/*_1.fq.gz
do
    SAMPLE=$(basename ${i} _1.fq.gz)

    echo "Trimming: " ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz

    /programs/bbmap-38.96/bbmerge.sh \
        in1=$FASTQ_TRIM/${SAMPLE}_1.fq.gz  \
        in2=$FASTQ_TRIM/${SAMPLE}_2.fq.gz \
        out=$MERGE_DIR/${SAMPLE}_trimmed_bbmerge.fq.gz \
        ihist=$OUT_DIR_INSERT_SIZE/${SAMPLE}_trimmed_bbmerge.txt
done



