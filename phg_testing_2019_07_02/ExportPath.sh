docker1 run --name cbsu_phg_container_exportPath --rm \
        -v /workdir/mbb262/PHGDocumentationTest/DockerOutput/:/tempFileDir/outputDir/ \
        -v /workdir/mbb262/PHGDocumentationTest/configSQLiteDocker.txt:/tempFileDir/data/configSQLiteDocker.txt \
        -t maizegenetics/phg \
        /ExportPath.sh configSQLiteDocker.txt CONSENSUS testOutput1.vcf