@file:Repository("http://maven.imagej.net/content/groups/public")
@file:DependsOn("com.github.samtools:htsjdk:2.19.0")
@file:DependsOn("net.maizegenetics:tassel:5.2.60")
@file:DependsOn("org.jetbrains.kotlinx:kotlinx-cli-jvm:0.3.5")

import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.ValidationStringency
import kotlinx.cli.*
import net.maizegenetics.util.Utils
import java.io.File
import kotlin.system.exitProcess

/**
 * data class that holds read counts for each parent and for ambiguous reads
 */
data class hybridReadMaps(val map1: Map<String, Int>, val map2: Map<String, Int>, val mapBoth: Map<String, Int>)

/**
 * Function reads in files and writes results table to output
 */
fun buildHybridCountMatrixFile(firstFile: String, secondFile: String, suffix1: String, suffix2: String, inputTranscriptListFile: String, outputFileName: String, minAlignLengthProp: Double, laxMaxNMProp: Double, strictMaxNMProp: Double, minMapQ: Int = 48) {
    // get transcript names from file
    val transcriptNames = Utils.getBufferedReader(inputTranscriptListFile).readLines().sorted()

    Utils.getBufferedWriter(outputFileName).use { output ->
        // add header
        output.write("transcript_id\tparent1_count\tparent2_count\tboth_equal_count\n")

        // make maps from files
        val maps = computeHybridReadCounts(File(firstFile), File(secondFile), suffix1, suffix2, minAlignLengthProp, laxMaxNMProp, strictMaxNMProp, minMapQ)

        //print counts to files
        transcriptNames.forEach {
            output.write("$it\t${maps.map1[it]?:0}\t${maps.map2[it]?:0}\t${maps.mapBoth[it]?:0}\n")
        }
    }

}

/**
 * Function counts reads and separates mapping by the best parent
 */
fun computeHybridReadCounts(firstFile: File, secondFile: File, suffix1: String, suffix2: String, minAlignLengthProp: Double, laxMaxNMProp: Double, strictMaxNMProp: Double, minMapQ: Int = 48): hybridReadMaps{
    // check that parameters are within bounds
    check(minAlignLengthProp in 0.0 .. 1.0) {"Error, minAlignmentLengthProportion is not between 0 and 1.0"}
    check(laxMaxNMProp in 0.0 .. 1.0) {"Error, laxMaxNMProp is not between 0 and 1.0"}
    check(strictMaxNMProp in 0.0 .. 1.0) {"Error, strictMaxNMProp is not between 0 and 1.0"}

    // make maps to store counts
    val transcriptMap1 = mutableMapOf<String, Int>()
    val transcriptMap2 = mutableMapOf<String, Int>()
    val transcriptMapBoth = mutableMapOf<String, Int>()

    // open up readers for both parent 1 and 2 sam files
    val reader1 = SamReaderFactory.makeDefault()
        .validationStringency(ValidationStringency.SILENT)
        .open(firstFile)

    val reader2 = SamReaderFactory.makeDefault()
        .validationStringency(ValidationStringency.SILENT)
        .open(secondFile)

    // and iterators
    val iterator1 = reader1.iterator()
    val iterator2 = reader2.iterator()

    var currentReadName = ""
    var currentRecords1 = mutableListOf<SAMRecord>()
    var currentRecords2 = mutableListOf<SAMRecord>()

    // large loop handles iteration through file 1
    while(iterator1.hasNext()) {
        val currentRecord1 = iterator1.next()

        // as long as we keep seeing records from the same read, add them to the records list
        if(currentRecord1.readName == currentReadName) {
            currentRecords1.add(currentRecord1)
        } else {
            // we have found all the records from parent1 with this read name
            // now we have to find all the records from parent 2

            // print statements for debugging
            //println("Found all records for sam1")
            //currentRecords1.forEach{ println(it)}

            // small loop iterates through file 2 until we "catch up" to file 1's location
            while(iterator2.hasNext()) {
                val currentRecord2 = iterator2.next()

                // add records from current read to the list
                if (currentRecord2.readName == currentReadName) {
                    currentRecords2.add(currentRecord2)
                } else {
                    // we have found all the records from parent2 with this read name
                    // now we have to filter the records

                    // print statements for debugging
                    //println("Found all records for sam2")
                    //currentRecords2.forEach{ println(it)}

                    // first round of filters is lax, so we catch subpar alignments to both parents
                    val passedAlignments1 = currentRecords1.filter { passesAlignmentFilter(it, minAlignLengthProp, laxMaxNMProp, minMapQ) }
                    val passedAlignments2 = currentRecords2.filter { passesAlignmentFilter(it, minAlignLengthProp, laxMaxNMProp, minMapQ) }

                    // if we have passing alignments to both parents
                    if(passedAlignments1.isNotEmpty() && passedAlignments2.isNotEmpty()) {

                        // based on minimap's ordering, the first alignments in the list should be the "best"
                        val best1 = passedAlignments1.first()
                        val best2 = passedAlignments2.first()

                        // decide on the better alignment, based on edit distance
                        if (best1.contig.removeSuffix(suffix1) == best2.contig.removeSuffix(suffix2)) {
                            val nm1 = best1.getIntegerAttribute("NM")
                            val nm2 = best2.getIntegerAttribute("NM")
                            if (nm1 < nm2) {
                                //println("sam1 has the better match")
                                if (passesAlignmentFilter(best1, minAlignLengthProp, strictMaxNMProp, minMapQ)) { transcriptMap1[best1.contig.removeSuffix(suffix1)] = (transcriptMap1[best1.contig.removeSuffix(suffix1)]?:0) + 1 }
                            } else if (nm1 > nm2) {
                                //println("sam2 has the better match")
                                if (passesAlignmentFilter(best2, minAlignLengthProp, strictMaxNMProp, minMapQ)) { transcriptMap2[best2.contig.removeSuffix(suffix2)] = (transcriptMap2[best2.contig.removeSuffix(suffix2)]?:0) + 1 }
                            } else { // nm1 == nm2
                                //println("both matches equal")
                                if (passesAlignmentFilter(best1, minAlignLengthProp, strictMaxNMProp, minMapQ) && passesAlignmentFilter(best2, minAlignLengthProp, strictMaxNMProp, minMapQ)) {
                                    transcriptMapBoth[best1.contig.removeSuffix(suffix1)] = (transcriptMapBoth[best1.contig.removeSuffix(suffix1)]?:0) + 1
                                    }
                            }
                        }  //else { println("best reads on different transcripts") }
                    } //else { println("No records passed the filter") }

                    // get the next read's name
                    currentReadName = currentRecord1.readName
                    //println("Next read: $currentReadName")

                    // this script only works if the two sam files have the same reads in the same order, and at least one record per read, including unmapped
                    if (currentRecord2.readName != currentReadName) {
                        println("ERROR: read $currentReadName not found in file 2. Check that read orders are consistent between sam files")
                        exitProcess(1)
                    }

                    currentRecords1 = mutableListOf(currentRecord1)
                    currentRecords2 = mutableListOf(currentRecord2)
                    //println()
                   break
                }
            }
        }
    }
    // catching the last read(s) in the file
    while(iterator2.hasNext()) {
        val currentRecord2 = iterator2.next()

        if (currentRecord2.readName == currentReadName) {
            currentRecords2.add(currentRecord2)
        } else {
            break
        }
    }
    // lax filter
    val passedAlignments1 = currentRecords1.filter { passesAlignmentFilter(it, minAlignLengthProp, laxMaxNMProp, minMapQ) }
    val passedAlignments2 = currentRecords2.filter { passesAlignmentFilter(it, minAlignLengthProp, laxMaxNMProp, minMapQ) }

    // get best alignment
    if(passedAlignments1.isNotEmpty() && passedAlignments2.isNotEmpty()) {
        val best1 = passedAlignments1.first()
        val best2 = passedAlignments2.first()

        if (best1.contig.removeSuffix(suffix1) == best2.contig.removeSuffix(suffix2)) {
            val nm1 = best1.getIntegerAttribute("NM")
            val nm2 = best2.getIntegerAttribute("NM")
            if (nm1 < nm2) {
                if (passesAlignmentFilter(best1, minAlignLengthProp, strictMaxNMProp, minMapQ)) { transcriptMap1[best1.contig.removeSuffix(suffix1)] = (transcriptMap1[best1.contig.removeSuffix(suffix1)]?:0) + 1 }
            } else if (nm1 > nm2) {
                if (passesAlignmentFilter(best2, minAlignLengthProp, strictMaxNMProp, minMapQ)) { transcriptMap2[best2.contig.removeSuffix(suffix2)] = (transcriptMap2[best2.contig.removeSuffix(suffix2)]?:0) + 1 }
            } else { // nm1 == nm2
                if (passesAlignmentFilter(best1, minAlignLengthProp, strictMaxNMProp, minMapQ) && passesAlignmentFilter(best2, minAlignLengthProp, strictMaxNMProp, minMapQ)) {
                    transcriptMapBoth[best1.contig.removeSuffix(suffix1)] = (transcriptMapBoth[best1.contig.removeSuffix(suffix1)]?:0) + 1
                }
            }
        }
    }
    return hybridReadMaps(transcriptMap1, transcriptMap2, transcriptMapBoth)
}

/**
 * Function to check to see if an alignment record will pass the requested filters.
 */
fun passesAlignmentFilter(currentRecord: SAMRecord, minAlignLengthProp: Double, maxNMProp: Double, minMapQ : Int = 48) : Boolean {
    val alignReadLength = currentRecord.cigar.cigarElements.filter { it.operator.consumesReadBases() }.filter { !it.operator.isClipping }.map { it.length }.sum()
    val actualReadLength = currentRecord.cigar.filter { when(it.operator) {
        CigarOperator.M, CigarOperator.I, CigarOperator.S, CigarOperator.EQ, CigarOperator.X, CigarOperator.H -> true
        else -> false
    } }.map { it.length }.sum()

    val editDist = currentRecord.getIntegerAttribute("NM")
    val mapQ = currentRecord.mappingQuality

    //println("is mapped ${!currentRecord.readUnmappedFlag}, is not clipped ${!currentRecord.cigar.isClipped()}, mapQ $mapQ, alignment length prop ${(alignReadLength.toDouble()/actualReadLength)}, nm prop ${(editDist.toDouble()/alignReadLength)}")

    //Checking both alignmentLengthProp and NM Prop
    return (!currentRecord.readUnmappedFlag && !currentRecord.cigar.isClipped() && mapQ > minMapQ  && (alignReadLength.toDouble()/actualReadLength) > minAlignLengthProp && (editDist.toDouble()/alignReadLength) < maxNMProp)
}

// define command line arguments
val parser = ArgParser("hybrid read counter")
val inFile1 by parser.option(ArgType.String, shortName = "i", description = "input file - sam aligned to parent 1").required()
val inFile2 by parser.option(ArgType.String, shortName = "j", description = "input file - sam aligned to parent 2").required()
val suffix1 by parser.option(ArgType.String, shortName = "p", description = "suffix appended to transcript names in sam of parent 1").default("")
val suffix2 by parser.option(ArgType.String, shortName = "q", description = "suffix appended to transcript names in sam of parent 2").default("")
val transcriptList by parser.option(ArgType.String, shortName = "t", description = "transcript list").required()
val outputFile by parser.option(ArgType.String, shortName = "o", description = "output table").required()
val minAlignLengthProp by parser.option(ArgType.Double, shortName = "a", description = "minimum proportion of read length aligned").default(0.9)
val laxMaxNMProp by parser.option(ArgType.Double, shortName = "l", description = "maximum proportion of read length NM lax boundary").default(0.1)
val strictMaxNMProp by parser.option(ArgType.Double, shortName = "m", description = "maximum proportion of read length NM strict boundary").default(0.02)
val minMapQ by parser.option(ArgType.Int, shortName = "n", description = "minimum mapQ score").default(48)
parser.parse(args)

buildHybridCountMatrixFile(inFile1, inFile2, suffix1, suffix2, transcriptList, outputFile, minAlignLengthProp, laxMaxNMProp, strictMaxNMProp, minMapQ)
