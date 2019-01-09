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

samples_ch.into { samples_ch1; samples_ch2 }




/* fastQC */


process rawFastQC {

	cpus = 4 // set as option --thread
	memory = '1GB' //as determined by fastqc - 250mb per core
	
	input:
    	set sampleid, uuid, file(fq1), file(fq2) from samples_ch1
	
	output:
		file "${uuid}_raw_fastqc_log.txt" 
		file "*fastqc.*"
		set uuid, file(fq1), file(fq2) into fastqc_out_ch
	
	tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/fastqc", mode: 'copy', pattern: "${uuid}*"
	
	"""
	cat $fq1 $fq2 > ${uuid}_raw.fq.gz
	fastqc --threads 4 ${uuid}_raw.fq.gz > ${uuid}_raw_fastqc_log.txt
	rm ${uuid}_raw.fq.gz
	"""

	
}


process bbDuk {

	cpus = 4 //see threads=4
	memory = '4GB' //see -Xmx4g
	
	input:
		set sampleid, uuid, file(fq1), file(fq2) from samples_ch2
	
	output:
		set uuid, file("${uuid}_clean.1.fq.gz"), file("${uuid}_clean.2.fq.gz") into bbduk_out_ch
	
	tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/clean_fastq", mode: 'copy'
	
	"""
	bbduk.sh in1=$fq1 in2=$fq2 out1=${uuid}_clean.1.fq out2=${uuid}_clean.2.fq \
				ref=$params.bbduk_adapaters ktrim=r k=23 mink=11 hdist=1 tpe tbo -Xmx4g threads=4
	gzip ${uuid}_clean.1.fq ${uuid}_clean.2.fq
	"""
}

bbduk_out_ch.into { bbduk_out_ch1; bbduk_out_ch2 }

/*
process spades {
	
	input:
		set uuid, file("${uuid}_clean.1.fq.gz"), file("${uuid}_clean.2.fq.gz") from bbduk_out_ch1
	
	output:
		file "${uuid}_*"
		
	tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/spades", mode: 'copy', pattern: "${uuid}_*"
	
	"""
	spades.py --careful -o spades -1 ${uuid}_clean.1.fq.gz -2 ${uuid}_clean.2.fq.gz
	cp spades/contigs.fasta ${uuid}_spades_contigs.fa
	cp spades/assembly_graph.fastg ${uuid}_spades_assembly_graph.fastg
	cp spades/assembly_graph_with_scaffolds.gfa ${uuid}_spades_assembly_graph_with_scaffolds.gfa
	cp spades/spades.log ${uuid}_spades.log
	
	"""
	
}
*/

process velvet {

	memory = '12GB'
	cpus = 1
	echo true

	input:
		set uuid, file("${uuid}_clean.1.fq.gz"), file("${uuid}_clean.2.fq.gz") from bbduk_out_ch2
	
	output:
		file "${uuid}_*"

		
	tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/velvet", mode: 'copy', pattern: "${uuid}_*"
	
	"""
	VelvetOptimiser.pl -s 53 -e 53 -x 4 -f '-shortPaired -fastq.gz -separate ${uuid}_clean.1.fq.gz ${uuid}_clean.2.fq.gz' -t 1
	cp auto_data_*/contigs.fa ${uuid}_velvet_contigs.fa
	cp auto_data_*/stats.txt ${uuid}_velvet_stats.txt
	cp *Logfile.txt ${uuid}_velvet_logfile.txt
	"""

}
