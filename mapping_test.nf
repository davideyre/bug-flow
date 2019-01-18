#!/usr/bin/env nextflow


/* test pipeline structure */

/* see example here - https://github.com/cbcrg/grape-nf/blob/master/grape.nf#L191 */

/*** TO DO
 - add commands to run unicycler+spades
 - generate guids and log of these during pipeline
 ***/

def getShortId( str ) {
  return str.substring(0,8) 
}


/* parameters */
params.index = "example_data/file_list.csv"
params.outputPath = "example_output"
params.refFile = "example_data/R00000419.fasta"
params.threads = 4

/* initial logging */
log.info "Pipeline Test -- version 0.1"
log.info "Input file              : ${params.index}"
log.info "Output path             : ${params.outputPath}"

/* input validation */
outputPath = file(params.outputPath)

Channel
    .fromPath(params.index)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleid, row.uuid, file(row.fq1), file(row.fq2)) }
    .set { samples_ch }

samples_ch.into { samples_ch1; samples_ch2 }

refFasta = file(params.refFile)
threads = params.threads


// Build BWA index for reference fasta file
process bwa_index {
  
    input:
        file refFasta
	
	output:
		file "*" into bwa_index
	
	
	tag {refFasta}
	publishDir "$outputPath/reference", mode: 'copy'

	
    """
    bwa index ${refFasta}
    """
}


// Map reads to reference genome with BWA MEM
process bwa{

    input:
    	set sampleid, uuid, file(fq1), file(fq2) from samples_ch1
    	file "*" from bwa_index
    	file refFasta

    output:
    	set uuid, file("${uuid}_alignment.sam") into bwa_mapped
    
    
    tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/bwa_mapped", mode: 'copy'

	// -R option specifies read group header line for output
	// 		chosen settings based on those used historically
    """
    bwa mem -R '@RG\tID:${uuid}\tSM:null\tLB:null\tCN:null' \
    		-t ${threads} \
    		${refFasta} \
    		${fq1} \
    		${fq2} \
    > ${uuid}_alignment.sam
    """
}
