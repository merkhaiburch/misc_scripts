/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/output_counts/goodman_random_count/goodman_debug_chrom1gshoot.txt \
    -Xmx200g \
    -CountAssociationsPlugin \
    -intervals /workdir/mbb262/genic_intergenic_intervals_b73_v4.49.csv \
    -gwasResults /workdir/mbb262/results/kremling/processed_results/chrom_1_fast_assoc_results_randomSNPs_GShoot_Kremling_2018_randomSNPs_reformatted_zack.csv \
    -outputCountFile /workdir/mbb262/output_counts/goodman_kremling_random_count/outputCountRandom_chrom_1_fast_assoc_results_randomSNPs_GShoot_Kremling_2018_randomSNPs_reformatted_zack.csv.txt \
    -endPlugin


/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/output_counts/goodman_random_count/goodman_debug.txt \
    -Xmx200g \
    -CountAssociationsPlugin \
    -intervals /workdir/mbb262/genic_intergenic_intervals_b73_v4.49.csv \
    -gwasResults /workdir/mbb262/results/kremling/processed_results/chrom_1_fast_assoc_results_randomSNPs_L3Base_Kremling_2018_randomSNPs_reformatted_zack.csv \
    -outputCountFile /workdir/mbb262/output_counts/goodman_random_count/outputCountRandom_chrom_1_fast_assoc_results_randomSNPs_L3Base_Kremling_2018_randomSNPs_reformatted_zack.csv.txt \
    -endPlugin

/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug /workdir/mbb262/output_counts/goodman_random_count/goodman_debug.txt \
    -Xmx200g \
    -CountAssociationsPlugin \
    -intervals /workdir/mbb262/genic_intergenic_intervals_b73_v4.49.csv \
    -gwasResults /workdir/mbb262/results/kremling/processed_results/chrom_5_fast_assoc_results_randomSNPs_LMAD_Kremling_2018_randomSNPs_reformatted_zack.csv \
    -outputCountFile /workdir/mbb262/output_counts/goodman_random_count/outputCountRandom_chrom_5_fast_assoc_results_randomSNPs_LMAD_Kremling_2018_randomSNPs_reformatted_zack.csv.txt \
    -endPlugin
