@file:DependsOn("com.github.samtools:htsjdk:2.19.0")
@file:DependsOn("net.maizegenetics:tassel:5.2.60")

import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.ValidationStringency
import net.maizegenetics.util.Utils
import java.io.File


/**
 * Method to take in a SAMDirectory and for each SAM file apply filters to each alignment and then sum up the number of reads which pass the filter for each transcript.
*/
fun buildCountMatrixFile(samDir : String, inputTranscriptListFile: String ,outputFileName:String ,minAlignLengthProp: Double, maxNMProp: Double, minMapQ: Int = 48) {
    val transcriptNames = Utils.getBufferedReader(inputTranscriptListFile).readLines().sorted()

    Utils.getBufferedWriter(outputFileName).use { output ->
            //add header
            output.write("\t${transcriptNames.joinToString("\t")}\n")

            val samFiles = File(samDir).walk()
                .filter { it.isFile }
                .filter { !it.isHidden }
                .map { Pair(it.nameWithoutExtension, computeReadCounts(it, minAlignLengthProp, maxNMProp, minMapQ)) }
                .forEach { readCounts ->
                    val taxaName = readCounts.first
                    output.write("${taxaName}\t${transcriptNames.map { readCounts.second[it]?:0 }.joinToString("\t")}\n")
                }
        }
}

/**
 * Function to count the number of reads which align to each transcript which pass the input filters.
*/
fun computeReadCounts(samFile: File, minAlignLengthProp: Double, maxNMProp: Double, minMapQ: Int = 48) : Map<String,Int>{
        check(minAlignLengthProp in 0.0 .. 1.0) {"Error, minAlignmentLengthProportion is not between 0 and 1.0"}
        check(maxNMProp in 0.0 .. 1.0) {"Error, maxNMProp is not between 0 and 1.0"}
        val transcriptMap = mutableMapOf<String, Int>()

        val reader = SamReaderFactory.makeDefault()
            .validationStringency(ValidationStringency.SILENT)
            .open(samFile)

        val samIterator = reader.iterator()

        var currentReadName = ""
        var currentRecords = mutableListOf<SAMRecord>()

        while(samIterator.hasNext()) {
            val currentRecord = samIterator.next()

            //Skip over if its unmapped
            if(currentRecord.readUnmappedFlag) {
                continue
            }

            if(currentRecord.readName != currentReadName) {
                val passedAlignments = currentRecords.filter { passesAlignmentFilter(it, minAlignLengthProp, maxNMProp, minMapQ) }

                if(passedAlignments.isNotEmpty()) {
                    transcriptMap[passedAlignments.first().contig] = (transcriptMap[passedAlignments.first().contig]?:0) + 1
                }

                //reset the readname
                currentReadName = currentRecord.readName
                currentRecords = mutableListOf()
            }
            currentRecords.add(currentRecord)

        }
        reader.close()

        return transcriptMap
}

/**
 * Function to check to see if an alignment record will pass the requested filters.
*/
fun passesAlignmentFilter(currentRecord: SAMRecord, minAlignLengthProp: Double, maxNMProp: Double, minMapQ : Int = 48) : Boolean {
        check(minAlignLengthProp in 0.0 .. 1.0) {"Error, minAlignmentLengthProportion is not between 0 and 1.0"}
        check(maxNMProp in 0.0 .. 1.0) {"Error, maxNMProp is not between 0 and 1.0"}

        val alignReadLength = currentRecord.cigar.cigarElements.filter { it.operator.consumesReadBases() }.filter { !it.operator.isClipping }.map { it.length }.sum()
        val actualReadLength = currentRecord.cigar.filter { when(it.operator) {
            CigarOperator.M, CigarOperator.I, CigarOperator.S, CigarOperator.EQ, CigarOperator.X, CigarOperator.H -> true
            else -> false
        } }.map { it.length }.sum()


        val editDist = currentRecord.getIntegerAttribute("NM")
        val mapQ = currentRecord.mappingQuality
        //Checking both alignmentLengthProp and NM Prop
        return (!currentRecord.cigar.isClipped() && mapQ > minMapQ  && (alignReadLength.toDouble()/actualReadLength) > minAlignLengthProp && (editDist.toDouble()/alignReadLength) < maxNMProp)
}


// Unit test
// val fullTranscriptSamDirUnit = "/home/jupyter-mbb262/data/merritt/sam_unit/"

// val transcriptListUnit = "/home/jupyter-mbb262/data/merritt/transcript_ids.txt"

// val outputFullTranscriptFileUnit = "/home/jupyter-mbb262/data/merritt/counts_output/unit_test_count.txt"

// 0.02 = only allow 3 mismatches, need a mapq of at least 49

// buildCountMatrixFile(fullTranscriptSamDirUnit, transcriptListUnit, outputFullTranscriptFileUnit, .9, 0.02, 47)



// Biological samples
val fullTranscriptSamDir = "/workdir/mbb262/nam_hybrid_rnaseq/output/minimap_alignments/"

val outputFullTranscriptFile = "/workdir/mbb262/nam_hybrid_rnaseq/output/counts/nam_inbred_hybrid_minimap_count.txt"

val transcriptList = "/workdir/mbb262/nam_hybrid_rnaseq/output/counts/all_nam_transcript_ids.txt"

buildCountMatrixFile(fullTranscriptSamDir, transcriptList, outputFullTranscriptFile, .9, 0.02, 48)
