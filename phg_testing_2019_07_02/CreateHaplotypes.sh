taxonList=(LineA1 LineA LineB1 LineB RefA1 Ref)

REF_DIR=/workdir/mbb262/PHGDocumentationTest/ref/
FASTQ_DIR=/workdir/mbb262/PHGDocumentationTest/WGSFastq
DB=/workdir/mbb262/PHGDocumentationTest/DockerOutput/phgSmallSeq.db
CONFIG_FILE=/workdir/mbb262/PHGDocumentationTest/configSQLiteDocker.txt
CONFIG_FILE_IN_DOCKER=/tempFileDir/data/configSQLiteDocker.txt
GVCF_OUTPUT_DIR=/workdir/mbb262/PHGDocumentationTest/gvcfOut

for TAXON in "${taxonList[@]}"
do
docker1 run --name small_seq_test_container --rm \
-w / \
-v ${REF_DIR}:/tempFileDir/data/reference/ \
-v ${FASTQ_DIR}:/tempFileDir/data/fastq/ \
-v ${DB}:/tempFileDir/outputDir/phgSmallSeq.db \
-v ${CONFIG_FILE}:${CONFIG_FILE_IN_DOCKER} \
-v ${GVCF_OUTPUT_DIR}:/tempFileDir/data/gvcfs/ \
-t maizegenetics/phg /CreateHaplotypesFromFastq.groovy -config ${CONFIG_FILE_IN_DOCKER} -t ${TAXON} -op single -fastqs ${TAXON}_R1.fastq

done


