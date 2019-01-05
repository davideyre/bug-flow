#!/usr/bin/env nextflow

/* test pipeline structure */

/* see example here - https://github.com/cbcrg/grape-nf/blob/master/grape.nf#L191 */

/*** TO DO
 - allow multiple input files to be read
 - add bbduk for adapter trimming (+/- quality trimming)
 - add kmergenie for velvet kmer finder
 - add commands to run vanilla velvet (with kmer size), spades, unicycler+spades
 - generate guids, log of these, and output file structure
 ***/

/* parameters */
params.dataPath = "example_data"
params.outputPath = "example_output"

/* initial logging */
log.info "Pipeline Test -- version 0.1"
log.info "Data path               : ${params.dataPath}"
log.info "Output path             : ${params.outputPath}"

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
		file "log.txt" 
		file "*fastqc.*" 
	
	"""
	fastqc $fastq_file  > log.txt
	"""

}

/*
process kmerGenie { 
	

	'''
	source activate py27
	kmergenie
	source deactivate
	'''

}

*/

