#!/usr/bin/env nextflow

/* test pipeline structure */

/* see example here - https://github.com/cbcrg/grape-nf/blob/master/grape.nf#L191 */

/*** TO DO
 - allow multiple input files to be read
 - add ps, date, sed, grep, egrep, awk to docker image
 ***/

/* parameters */
params.dataPath = "example_data"
params.outputPath = "example_output"

/* initial logging */
log.info "Pipeline Test -- version 0.1"
log.info "Data path               : ${params.dataPath}"
log.info "Output path               : ${params.outputPath}"

/* input validation */
dataPath = file(params.dataPath)
outputPath = file(params.outputPath)
fastq_file_location = file("$dataPath/IR72_1.fastq.gz")

/* fastQC */
process fastQC {

	publishDir outputPath, mode: 'copy'
	
	input:
	    file fastq_file from fastq_file_location
	
	output:
		file "log.txt" into fastQC_stdout
		file "*fastqc.*" into fastQC_outputFiles
	
	"""
	fastqc $fastq_file  > log.txt
	"""

}

