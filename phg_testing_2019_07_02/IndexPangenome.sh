DB=/workdir/mbb262/PHGDocumentationTest/DockerOutput/phgSmallSeq.db
OUTPUT_DIR=/workdir/mbb262/PHGDocumentationTest/DockerOutput/
CONFIG_FILE=/workdir/mbb262/PHGDocumentationTest/configSQLiteDocker.txt
docker1 run --name index_pangenome_container --rm \
                    -w / \
                    -v ${OUTPUT_DIR}:/tempFileDir/outputDir/pangenome/ \
                    -v ${DB}:/tempFileDir/outputDir/phgSmallSeq.db \
                    -v ${CONFIG_FILE}:/tempFileDir/data/configSQLiteDocker.txt \
                    -t maizegenetics/phg \
                    /IndexPangenome.sh pangenome_fasta configSQLiteDocker.txt CONSENSUS 4G 15 10