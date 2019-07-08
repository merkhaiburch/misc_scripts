docker1 run --name phg_container_consensus --rm \
        -v /workdir/mbb262/PHGDocumentationTest/ref/:/tempFileDir/data/reference/ \
        -v /workdir/mbb262/PHGDocumentationTest/DockerOutput/:/tempFileDir/outputDir/ \
        -v /workdir/mbb262/PHGDocumentationTest/configSQLiteDocker.txt:/tempFileDir/data/configSQLiteDocker.txt \
        -t maizegenetics/phg \
        /CreateConsensi.sh /tempFileDir/data/configSQLiteDocker.txt Ref.fa GATK_PIPELINE CONSENSUS