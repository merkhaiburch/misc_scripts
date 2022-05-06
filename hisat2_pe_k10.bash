


# Map all reads to b73xky21 --> k10
for i in $FASTQC_TRIM_DIR/*_1.fq.gz
do
    SAMPLE=$(basename ${i} _1.fq.gz)

    echo "Aligning: " ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz

    hisat2 \
        -p $N_THREADS \
        -x $PROJ_DIR/hisat2_index \
        -k 10 \
        --no-softclip \
        -1 $FASTQC_TRIM_DIR/${SAMPLE}_1.fq.gz \
        -2 $FASTQC_TRIM_DIR/${SAMPLE}_2.fq.gz \
        -S $OUT_PATH/${SAMPLE}_b73ky21_k10_output.sam --new-summary \
        2> $PROJ_DIR/output/summaries/${SAMPLE}_b732ky21_k10_summary.txt
done


# Map all reads to b73 --> k10
for i in $FASTQC_TRIM_DIR/*_1.fq.gz
do
    SAMPLE=$(basename ${i} _1.fq.gz)

    echo "Aligning: " ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz

    hisat2 \
        -p $N_THREADS \
        -x $PROJ_DIR/hisat2_index_b73 \
        -k 10 \
        --no-softclip \
        -1 $FASTQC_TRIM_DIR/${SAMPLE}_1.fq.gz \
        -2 $FASTQC_TRIM_DIR/${SAMPLE}_2.fq.gz \
        -S $OUT_PATH/${SAMPLE}_b73_k10_output.sam --new-summary \
        2> $PROJ_DIR/output/summaries/${SAMPLE}_b73_k10_summary.txt
done

# Map all reads to ky21 --> k10
for i in $FASTQC_TRIM_DIR/*_1.fq.gz
do
    SAMPLE=$(basename ${i} _1.fq.gz)

    echo "Aligning: " ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz

    hisat2 \
        -p $N_THREADS \
        -x $PROJ_DIR/hisat2_index_ky21 \
        -k 10 \
        --no-softclip \
        -1 $FASTQC_TRIM_DIR/${SAMPLE}_1.fq.gz \
        -2 $FASTQC_TRIM_DIR/${SAMPLE}_2.fq.gz \
        -S $OUT_PATH/${SAMPLE}_ky21_k10_output.sam --new-summary \
        2> $PROJ_DIR/output/summaries/${SAMPLE}_ky21_k10_summary.txt
done