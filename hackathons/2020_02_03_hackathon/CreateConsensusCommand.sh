perl /workdir/mbb262/tassel-5-standalone/run_pipeline.pl -Xmx200g -debug -HaplotypeGraphBuilderPlugin -configFile /workdir/ag2484/ConfigFile -methods mummer4 -includeVariantContexts true -endPlugin \
                                                -RunHapConsensusPipelinePlugin -ref /workdir/ag2484/Zm-B73-REFERENCE-NAM-5.0.fa \
                                                                                -dbConfigFile /workdir/ag2484/ConfigFile \
                                                                                -collapseMethod CONSENSUS_mxdiv_0.001 \
                                                                                -collapseMethodDetails "\"CONSENSUS for creating Consensus\"" \
                                                                                -minFreq 0.01 \
                                                                                -endPlugin
