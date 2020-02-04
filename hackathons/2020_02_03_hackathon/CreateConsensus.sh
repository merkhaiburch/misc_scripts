dbConfigFile=

perl /workdir/mbb262/tassel-5-standalone/run_pipeline.pl $fullXmx \
  -debug \
  -HaplotypeGraphBuilderPlugin \
  -configFile $dbConfigFile \
  -methods ${HAPLOTYPE_METHOD} \
  -includeVariantContexts ${includeVariants} \
  -endPlugin \
  -RunHapConsensusPipelinePlugin \
  -ref ${REFERENCE_FILE} \
  -dbConfigFile ${dbConfigFile} \
  -collapseMethod ${COLLAPSE_METHOD} \
  -collapseMethodDetails "\"${COLLAPSE_METHOD} for creating Consensus\"" \
  -minFreq ${minFreq} \
  -endPlugin