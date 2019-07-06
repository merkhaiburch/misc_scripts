# This script assumes the config.txt file lives in the directory mounted to /tempFileDir/data
# It assumes the reference fasta (in this case, Zea_mays.AGPv4.dna.toplevelMtPtv3.fa) lives
# in the directory mounted to /tempFileDir/data/reference.
# It assumes the genome intervals file (in the example below, maizeRefAnchor_intervals_bed) lives
# in the directory mounted to /tempFileDir/answer.
# It assumes the genomeData file describing the reference (in the example below, B&3Ref_load_data.txt) lives
# in the directory mounted to /tempFileDir/data
# It assumes your sqlite database lives in the directory mounted to /tempFileDir/outputDir/  This in only relevant when running an SQLite database.  This path shows up in the config file, parameter "db".

# You must change "/workdir/user/DockerTuningTests/..." to match your own directory paths
docker1 run --name load_phg_container --rm \
        -v /workdir/mbb262/PHGDocumentationTest/DockerOutput/:/tempFileDir/outputDir/ \
        -v /workdir/mbb262/PHGDocumentationTest/ref/:/tempFileDir/data/reference/ \
        -v /workdir/mbb262/PHGDocumentationTest/:/tempFileDir/data/ \
        -v /workdir/mbb262/PHGDocumentationTest/:/tempFileDir/answer/ \
        -t maizegenetics/phg:latest \
        /LoadGenomeIntervals.sh configSQLiteDocker.txt Ref.fa intervals.bed load_data.txt true