#!/usr/bin/env nextflow

/* simple example pipeline for mapping and variant calling */

def getShortId( str ) {
  return str.substring(0,8) 
}


/* parameters */
params.index = "example_data/file_list.csv"
params.outputPath = "example_output"
params.refFile = "example_data/R00000419.fasta"
params.threads = 6

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
    samtools faidx ${refFasta}
    """
}


// Map reads to reference genome with BWA MEM
process bwa{

    input:
    	set sampleid, uuid, file(fq1), file(fq2) from samples_ch1
    	file "*" from bwa_index
    	file refFasta

    output:
    	set uuid, file("${uuid}.aligned.sam") into bwa_mapped
    
    
    tag "${getShortId(uuid)}"
	//publishDir "$outputPath/$uuid/bwa_mapped", mode: 'copy'

    """
    bwa mem -t ${threads} \
    		${refFasta} \
    		${fq1} \
    		${fq2} \
    > ${uuid}.aligned.sam
    """
}

//remove duplicates using samtools v 1.9
process removeDuplicates{

    input:
    	set uuid, file("${uuid}.aligned.sam") from bwa_mapped
    	file "*" from bwa_index

    output:
    	set uuid, file("${uuid}.bam"), file("${uuid}.bam.bai") into dup_removed
    
    tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/bwa_mapped", mode: 'copy', pattern: "${uuid}.ba*"

	//sort by name to run fixmate (to remove BWA artefacts) and provide info for markdup
	//sort by position to run markdup (and then remove duplicates)
    """
    samtools sort -@${threads} -n -o sorted.bam ${uuid}.aligned.sam
    samtools fixmate -m sorted.bam fixed.sorted.bam
    samtools sort -@${threads} -o fixed.resorted.bam fixed.sorted.bam
    samtools markdup -r fixed.resorted.bam ${uuid}.bam
    samtools index ${uuid}.bam
    """
}


//call SNPs using samtools
process snpCall{


    input:
    	set uuid, file("${uuid}.bam"), file("${uuid}.bam.bai") from dup_removed
    	file "*" from bwa_index
    	file refFasta
 
    output:
    	set uuid, file("${uuid}.vcf.gz") into snps_called
   
    tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/bwa_mapped", mode: 'copy'

	//use bcftools mpileup to generate vcf file
	//mpileup genearates the likelihood of each base at each site
	//call converts this to actual variants in the VCF file
	//filter applies filters - commented out for now
    """
    bcftools mpileup -f $refFasta ${uuid}.bam | \
    	bcftools call -Oz -mv --ploidy 1 --threads $threads > ${uuid}.vcf.gz
    """
    /*
        """
    bcftools mpileup -f $refFasta ${uuid}.bam | \
    	bcftools call -Ou -mv | \
    	bcftools filter -s LowQual -e '%QUAL<20 || DP>100' > ${uuid}.vcf
    """
    */

}


//filter SNPS


//generate concensus fasta file
process concensusFa{

	input:
		set uuid, file("${uuid}.vcf.gz") from snps_called
		file refFasta
	
	output:
		set uuid, file("${uuid}.fa") into fa_file
	
	tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/bwa_mapped", mode: 'copy'

	"""
	cat $refFasta | bcftools consensus ${uuid}.vcf.gz > ${uuid}.fa
	"""

}


