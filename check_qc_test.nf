#!/usr/bin/env nextflow

/* test pipeline structure */

/* see example here - https://github.com/cbcrg/grape-nf/blob/master/grape.nf#L191 */

/*** TO DO
 - allow multiple input files to be read
 - add ps, date, sed, grep, egrep, awk to docker image
 ***/

/* parameters */
params.dataPath = "/Users/davideyre/Drive/academic/infrastructure/nextflow/pipeline-test/example_data"

/* initial logging */
log.info "Pipeline Test -- version 0.1"
log.info "Data path               : ${params.dataPath}"

/* input validation */
dataPath = file(params.dataPath)
fastq_file_location = file("$dataPath/IR72_1.fastq.gz")

/* fastQC */
process fastQC {
	
	input:
	    file fastq_file from fastq_file_location
	
	output:
		file "log.txt" into fastQC_stdout
		file "*fastqc.*" into fastQC_outputFiles
	
	"""
	fastqc $fastq_file  > log.txt
	"""

}

/* save stdout from fastQC */
fastQC_stdout.subscribe { 
	it.copyTo(dataPath)
}

/* copy fastQC output files back to main directory */
fastQC_outputFiles
	.flatMap()
	.subscribe { it.copyTo(dataPath) }