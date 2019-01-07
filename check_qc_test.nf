#!/usr/bin/env nextflow


/* test pipeline structure */

/* see example here - https://github.com/cbcrg/grape-nf/blob/master/grape.nf#L191 */

/*** TO DO
 - add bbduk for adapter trimming (+/- quality trimming)
 - add kmergenie for velvet kmer finder
 - add commands to run vanilla velvet (with kmer size), spades, unicycler+spades
 - generate guids and log of these during pipeline
 ***/

def getShortId( str ) {
  return str.substring(0,8) 
}


/* parameters */
params.index = "example_data/file_list.csv"
params.outputPath = "example_output"

params.bbduk_adapaters = "/opt/conda/opt/bbmap-38.22-0/resources/adapters.fa"

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




/* fastQC */

process fastQC {

	//errorStrategy 'ignore' //set globally
	
	
	publishDir "$outputPath/$uuid/fastqc", mode: 'copy'
	
	input:
    	set sampleid, uuid, file(fq1), file(fq2) from samples_ch
	
	output:
		file "${uuid}_fastqc_log.txt" 
		file "*fastqc.*"
		set uuid, file(fq1), file(fq2) into fastqc_out_ch
	
	//tag each process with the guid
	tag "${getShortId(uuid)}"
	
	"""
	cat $fq1 $fq2 > ${uuid}_merge.fq.gz
	fastqc ${uuid}_merge.fq.gz > ${uuid}_fastqc_log.txt
	rm ${uuid}_merge.fq.gz
	"""

	
}

process bbDuk {
	
	publishDir "$outputPath/$uuid/clean_fastq", mode: 'copy'
	
	input:
		set uuid, file(fq1), file(fq2) from fastqc_out_ch
	
	output:
		set uuid, file("${uuid}_clean.1.fq.gz"), file("${uuid}_clean.2.fq.gz") into bbduk_out_ch
	
	tag "${getShortId(uuid)}"
	
	"""
	bbduk.sh in1=$fq1 in2=$fq2 out1=${uuid}_clean.1.fq out2=${uuid}_clean.2.fq \
				ref=$params.bbduk_adapaters ktrim=r k=23 mink=11 hdist=1 tpe tbo
	gzip ${uuid}_clean.1.fq ${uuid}_clean.2.fq
	"""
}

process kmerGenie {
	
	publishDir "$outputPath/$uuid/kmer_genie", mode: 'copy'
	
	echo true
	
	tag "${getShortId(uuid)}"
	
	input:
		set uuid, file("${uuid}_clean.1.fq.gz"), file("${uuid}_clean.2.fq.gz") from bbduk_out_ch
	
	output:
		file "kmergenie.log"
	
	"""
	ls -1 *.fq.gz > list_files
    kmergenie list_files -k 191 > kmergenie.log
	"""
	
}

