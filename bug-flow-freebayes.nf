#!/usr/bin/env nextflow


/* 
An example pipeline for mapping followed by variant calling and de novo assembly

Mapping - 
 - based on BWA mem
 - samtools / bcftools (NB. Indels are saved but not used for the consensus sequence)
 - filters
 
Assembly
 - spades 

*/

//function for coverting UUID to first 8 digits
def getShortId( str ) {
  if (input_type_is_csv) {
    return str.substring(0,8)
  }
  else {
    return str
  }
}


// parameters 
params.seqlist = "example_data/file_list.csv"
params.outputPath = "example_output"
params.refFile = "example_data/R00000419.fasta"
params.input_type_is_csv = true
params.runSpades = 1

// parameters if input_type_is_csv is false
params.indir = ""
params.readpat = ""

// initial logging
log.info "\n" 
log.info "BUGflow -- version 0.1"
log.info "Input sequence list    :  ${params.seqlist}"
log.info "Input reference file   :  ${params.refFile}"
log.info "Output path            :  ${params.outputPath}"
log.info "Container engine       :  ${workflow.containerEngine}"
log.info "Container              :  ${workflow.container}"
log.info "Profile                :  ${workflow.profile}"
log.info "\n"


// rename input parameters
refFasta = file(params.refFile)
outputPath = file(params.outputPath)
runSpades = params.runSpades

indir = params.indir
readpat = params.readpat
input_type_is_csv = params.input_type_is_csv
data_path = indir + readpat

//location for bbduk adapter sequences
bbduk_adapaters = "/opt/conda/opt/bbmap-38.22-0/resources/adapters.fa" //path within docker/singularity image

if (input_type_is_csv) {
  // set up initial channel based on CSV file
  Channel
      .fromPath(params.seqlist)
      .splitCsv(header:true)
      .map{ row-> tuple(row.uuid, file(row.fq1), file(row.fq2)) }
      .set { samples_ch }
}
else {
  // alternatively, set up channel from input_dir and readpat
  Channel
      .fromFilePairs(data_path, flat: true)
      .ifEmpty{ error "Cannot find any reads matching: ${params.readpat}" }
      .set { samples_ch }
}
samples_ch.into { samples_ch1; samples_ch2 }


// Build indexes for reference fasta file - bwa, samtools, repeatmask
process indexReference {
  
    input:
        file refFasta
	
	output:
		file "*" into ref_index
	
	tag {refFasta}
	publishDir "$outputPath/reference/${refFasta.baseName}", mode: 'copy'

    """
    #bwa index
    bwa index $refFasta
    
    #samtools index
    samtools faidx $refFasta
    
    #blast indexes for self-self blast
	makeblastdb -dbtype nucl -in $refFasta
    
    #reference mask
    genRefMask.py -r $refFasta -m 200 -p 95
    bgzip -c ${refFasta}.rpt.regions > ${refFasta.baseName}.rpt_mask.gz
	echo '##INFO=<ID=RPT,Number=1,Type=Integer,Description="Flag for variant in repetitive region">' > ${refFasta.baseName}.rpt_mask.hdr
	tabix -s1 -b2 -e3 ${refFasta.baseName}.rpt_mask.gz
    """
}



// initial fastQC
process rawFastQC {
	
	input:
        set uuid, file(fq1), file(fq2) from samples_ch1
	
	output:
		file "*"
	
	tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/raw_fastqc", mode: 'copy', pattern: "${uuid}*"
	
	"""
	cat $fq1 $fq2 > ${uuid}_raw.fq.gz
	fastqc --threads ${task.cpus} ${uuid}_raw.fq.gz > ${uuid}_raw_fastqc_log.txt
	rm ${uuid}_raw.fq.gz
	"""

}

//adapter trimming with bbDuk
process bbDuk {
	
	input:
		set uuid, file(fq1), file(fq2) from samples_ch2
	
	output:
		set uuid, file("${uuid}_clean.1.fq.gz"), file("${uuid}_clean.2.fq.gz") into bbduk_out_ch
	
	tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/clean_fastq", mode: 'copy'
	
	"""
	bbduk.sh in1=$fq1 in2=$fq2 out1=${uuid}_clean.1.fq out2=${uuid}_clean.2.fq \
				ref=$bbduk_adapaters ktrim=r k=23 mink=11 hdist=1 \
				tpe tbo -Xmx${task.memory.toGiga()}g threads=${task.cpus}
	gzip ${uuid}_clean.1.fq ${uuid}_clean.2.fq
	"""
}

//split cleaned reads into 3 channels - for repeat QC, assembly and mapping
bbduk_out_ch.into { bbduk_out_ch1; bbduk_out_ch2; bbduk_out_ch3 }


// repeat fastQC
process cleanFastQC {
	
	input:
    	set uuid, file("${uuid}_clean.1.fq.gz"), file("${uuid}_clean.2.fq.gz") from bbduk_out_ch1
	
	output:
		file "*"
	
	tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/clean_fastqc", mode: 'copy', pattern: "${uuid}*"
	
	"""
	cat ${uuid}_clean.1.fq.gz ${uuid}_clean.2.fq.gz > ${uuid}_clean.fq.gz
	fastqc --threads ${task.cpus} ${uuid}_clean.fq.gz > ${uuid}_clean_fastqc_log.txt
	rm ${uuid}_clean.fq.gz
	"""

}

if (runSpades == 1) {
	process spades {
	
	input:
		set uuid, file("${uuid}_clean.1.fq.gz"), file("${uuid}_clean.2.fq.gz") from bbduk_out_ch2
	
	output:
		file "${uuid}_*"
		
	tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/spades", mode: 'copy', pattern: "${uuid}_*"
	
	"""
	spades.py --careful -o spades -1 ${uuid}_clean.1.fq.gz -2 ${uuid}_clean.2.fq.gz \
		-t ${task.cpus} -m ${task.memory.toGiga()}
	cp spades/contigs.fasta ${uuid}_spades_contigs.fa
	cp spades/assembly_graph.fastg ${uuid}_spades_assembly_graph.fastg
	cp spades/assembly_graph_with_scaffolds.gfa ${uuid}_spades_assembly_graph_with_scaffolds.gfa
	cp spades/spades.log ${uuid}_spades.log
	"""
	}
}




// Map reads to reference genome with BWA MEM
process bwa {

    input:
    	set uuid, file("${uuid}_clean.1.fq.gz"), file("${uuid}_clean.2.fq.gz") from bbduk_out_ch3
    	file "*" from ref_index
    	file refFasta

    output:
    	set uuid, file("${uuid}.aligned.sam") into bwa_mapped
    
    tag "${getShortId(uuid)}"

	//don't add read group header here results in poorly formatted header
    """
    bwa mem -r 1.5 -O 6 \
    		-t ${task.cpus} \
    		$refFasta \
    		${uuid}_clean.1.fq.gz \
    		${uuid}_clean.2.fq.gz \
    > ${uuid}.aligned.sam
    """
}

//remove duplicates using samtools v 1.9
process removeDuplicates {

    input:
    	set uuid, file("${uuid}.aligned.sam") from bwa_mapped
    	file "*" from ref_index

    output:
    	set uuid, file("${uuid}.bam"), file("${uuid}.bam.bai") into dup_removed

    tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/bam", mode: 'copy', pattern: "${uuid}.ba*"

	//sort by name to run fixmate (to remove BWA artefacts) and provide info for markdup
	//sort by position to run markdup (and then remove duplicates)
    """
    samtools sort -@${task.cpus} -n -o sorted.bam ${uuid}.aligned.sam
    samtools fixmate -m sorted.bam fixed.sorted.bam
    samtools sort -@${task.cpus} -o fixed.resorted.bam fixed.sorted.bam
    samtools markdup -r fixed.resorted.bam ${uuid}.bam
    samtools index ${uuid}.bam
    """
}

process freebayes {
    input:
    	set uuid, file("${uuid}.bam"), file("${uuid}.bam.bai") from dup_removed
    	file refFasta
    	file "*" from ref_index
 
    output:
    	set uuid, file("raw_var.vcf") into raw_vcf
   
    tag "${getShortId(uuid)}"
    publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/freebayes", mode: 'copy'

	//generate a permissive list of all possible variants
 	"""
   	freebayes -p 1 --fasta-reference $refFasta ${uuid}.bam > raw_var.vcf
    """

}

//call SNPs using samtools call from mpileup file
process normalise_snps {

    input:
    	set uuid, file("raw_var.vcf") from raw_vcf
    	file refFasta
    	file "*" from ref_index
 
    output:
    	set uuid, file("${uuid}.bcf") into snps_called
   
    tag "${getShortId(uuid)}"
    publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/freebayes", mode: 'copy'
	
	
    """   
    bcftools norm -f $refFasta -m +any -Ou -o ${uuid}.bcf raw_var.vcf
    """

}


//filter SNPS
process filterSnps {

    input:
    	set uuid, file("${uuid}.bcf") from snps_called
    	file refFasta
    	file "*" from ref_index
 
    output:
    	set uuid, file("${uuid}.snps.vcf.gz"), file("${uuid}.snps.vcf.gz.csi") into filtered_snps
    	file "*"
   
    tag "${getShortId(uuid)}"
	//publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.snps.*"
	//publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.indels.*"
	//publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.zero_coverage.*"
	//publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.all.*"
	publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/freebayes", mode: 'copy'
	
	//use bcftools to filter normalised bcf file of variants from pileup and call
	//use one line for each filter condition and label
	//create index at end for random access and consensus calling
	
	//filters
		// quality >30
		// one read in each direction to support variant
		// not in a repeat region
		// consensus of >90% reads to support alternative allele
		// mask SNPs within 7 bp of INDEL
		// require high quality depth of 5 for call
	
    """
    #annotate vcf file with repetitive regions
	bcftools annotate -a ${refFasta.baseName}.rpt_mask.gz -c CHROM,FROM,TO,RPT \
		-h ${refFasta.baseName}.rpt_mask.hdr ${uuid}.bcf -Ob -o ${uuid}.masked.bcf.gz
    
    #filter vcf
    bcftools filter -S . -s Q30 -e '%QUAL<30' -Ou ${uuid}.masked.bcf.gz | \
    	bcftools filter -S . -s OneEachWay -e 'SAF == 0 || SAR ==0' -m+ -Ou | \
    	bcftools filter -S . -s RptRegion -e 'RPT=1' -m+ -Ou | \
    	bcftools filter -S . -s Consensus90 -e '((SAF+SAR)/(SRF+SRR+SAF+SAR))<=0.9' -m+ -Ou | \
    	bcftools filter -s SnpGap --SnpGap 7 -m+ -Ou | \
    	bcftools filter -S . -s SnpGap -e 'FILTER ~ "SnpGap"' -m+ -Ou | \
    	bcftools filter -S . -s HQDepth5 -e '(SAF+SAR)<=5' -m+ -Oz -o ${uuid}.all.vcf.gz
    
    
    #create vcf file with just SNPs
    bcftools filter -i 'TYPE="snp"' -m+ -Oz -o ${uuid}.snps.vcf.gz ${uuid}.all.vcf.gz 
    bcftools index ${uuid}.snps.vcf.gz
    
    #create vcf file with just INDELs
    bcftools filter -i 'TYPE!="snp"' -m+ -Oz -o ${uuid}.indels.vcf.gz ${uuid}.all.vcf.gz 
    bcftools index ${uuid}.indels.vcf.gz
    """


}

//generate consensus fasta file
process consensusFa {

	input:
		set uuid, file("${uuid}.snps.vcf.gz"), file("${uuid}.snps.vcf.gz.csi") from filtered_snps
		file refFasta
	
	output:
		set uuid, file("${uuid}.fa") into fa_file
	
	tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/fasta", mode: 'copy', pattern: "${uuid}.*"

	// call consensus sequence
		// -S flag in bcftools filter sets GT (genotype) to missing, with -M flag here
		// setting value to N
	"""
	cat $refFasta | bcftools consensus -H 1 -M "N" ${uuid}.snps.vcf.gz > ${uuid}.fa
	"""
}

