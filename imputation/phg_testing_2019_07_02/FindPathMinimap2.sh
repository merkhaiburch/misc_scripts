DB=/workdir/mbb262/PHGDocumentationTest/DockerOutput/phgSmallSeq.db
CONFIG_FILE=/workdir/mbb262/PHGDocumentationTest/configSQLiteDocker.txt
PANGENOME_DIR=/workdir/mbb262/PHGDocumentationTest/DockerOutput/

docker1 run --name small_seq_test_container --rm \
                    -w / \
                    -v /workdir/mbb262/PHGDocumentationTest/DockerOutput/:/tempFileDir/outputDir/pangenome/ \
                    -v /workdir/mbb262/PHGDocumentationTest/WGSFastq/:/tempFileDir/data/fastq/ \
                    -v /workdir/mbb262/PHGDocumentationTest/configSQLiteDocker.txt:/tempFileDir/data/configSQLiteDocker.txt \
                    -v ${DB}:/tempFileDir/outputDir/phgSmallSeq.db \
                    -v ${PANGENOME_DIR}:/tempFileDir/outputDir/pangenome/ \
                    -t maizegenetics/phg \
                    /FindPathMinimap2.sh pangenome_fasta configSQLiteDocker.txt \
                    CONSENSUS CONSENSUS,refRegionGroup \
                    HAP_COUNT_METHOD PATH_METHOD false
