#!/usr/bin/env nextflow

/* TO DO

set up memory and CPU for spades


*/


/* 
An example pipeline for mapping followed by variant calling and de novo assembly

Mapping - 
 - based on BWA mem
 - samtools / bcftools (NB. Indels are saved but not used for the consensus sequence)
 - filters
 
Assembly
 - spades or velvet optimiser based velvet

*/

//function for coverting UUID to first 8 digits
def getShortId( str ) {
  return str.substring(0,8) 
}


// parameters 
params.index = "example_data/file_list.csv"
params.outputPath = "example_output"
params.refFile = "example_data/R00000419.fasta"
params.threads = 6 //number of threads for multithreaded tasks
params.bbduk_adapaters = "/opt/conda/opt/bbmap-38.22-0/resources/adapters.fa" //path within docker image

// initial logging 
log.info "BUGflow -- version 0.1"
log.info "Input file              : ${params.index}"
log.info "Output path             : ${params.outputPath}"


// rename input parameters
refFasta = file(params.refFile)
threads = params.threads
outputPath = file(params.outputPath)


// set up initial channel based on CSV file
Channel
    .fromPath(params.index)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleid, row.uuid, file(row.fq1), file(row.fq2)) }
    .set { samples_ch }

samples_ch.into { samples_ch1; samples_ch2 }


// Build indexes for reference fasta file - bwa, samtools, repeatmask
process indexReference {
  
    input:
        file refFasta
	
	output:
		file "*" into ref_index
	
	cpus 1	
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

	cpus = threads // set as option --thread
	memory = "${threads*0.25} G" //as determined by fastqc - 250mb per core
	
	input:
    	set sampleid, uuid, file(fq1), file(fq2) from samples_ch1
	
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

	cpus = threads //see threads=4
	memory = 4.GB //see -Xmx4g, based on example runs
	
	input:
		set sampleid, uuid, file(fq1), file(fq2) from samples_ch2
	
	output:
		set uuid, file("${uuid}_clean.1.fq.gz"), file("${uuid}_clean.2.fq.gz") into bbduk_out_ch
	
	tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/clean_fastq", mode: 'copy'
	
	"""
	bbduk.sh in1=$fq1 in2=$fq2 out1=${uuid}_clean.1.fq out2=${uuid}_clean.2.fq \
				ref=$params.bbduk_adapaters ktrim=r k=23 mink=11 hdist=1 \
				tpe tbo -Xmx${task.memory.toGiga()}g threads=${task.cpus}
	gzip ${uuid}_clean.1.fq ${uuid}_clean.2.fq
	"""
}

//split cleaned reads into 3 channels - for repeat QC, assembly and mapping
bbduk_out_ch.into { bbduk_out_ch1; bbduk_out_ch2; bbduk_out_ch3 }


// repeat fastQC
process cleanFastQC {

	cpus = threads // set as option --thread
	memory = "${threads*0.25} G" //as determined by fastqc - 250mb per core
	
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

process spades {
	
	input:
		set uuid, file("${uuid}_clean.1.fq.gz"), file("${uuid}_clean.2.fq.gz") from bbduk_out_ch2
	
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


// process velvet {
// 
// 	memory = '12GB' //best guess from previous runs
// 	cpus = 1 //just use one for now given memory limitation
// 	echo true
// 
// 	input:
// 		set uuid, file("${uuid}_clean.1.fq.gz"), file("${uuid}_clean.2.fq.gz") from bbduk_out_ch2
// 	
// 	output:
// 		file "${uuid}_*"
// 
// 		
// 	tag "${getShortId(uuid)}"
// 	publishDir "$outputPath/$uuid/velvet", mode: 'copy', pattern: "${uuid}_*"
// 	
// 	"""
// 	VelvetOptimiser.pl -s 33 -e 171 -x 4 -f '-shortPaired -fastq.gz \
//		-separate ${uuid}_clean.1.fq.gz ${uuid}_clean.2.fq.gz' -t ${task.cpus}
// 	cp auto_data_*/contigs.fa ${uuid}_velvet_contigs.fa
// 	cp auto_data_*/stats.txt ${uuid}_velvet_stats.txt
// 	cp *Logfile.txt ${uuid}_velvet_logfile.txt
// 	"""
// 
// }


// Map reads to reference genome with BWA MEM
process bwa{

    input:
    	set uuid, file("${uuid}_clean.1.fq.gz"), file("${uuid}_clean.2.fq.gz") from bbduk_out_ch3
    	file "*" from ref_index
    	file refFasta

    output:
    	set uuid, file("${uuid}.aligned.sam") into bwa_mapped
    
    cpus threads
    tag "${getShortId(uuid)}"

    """
    bwa mem -t ${task.cpus} \
    		$refFasta \
    		${uuid}_clean.1.fq.gz \
    		${uuid}_clean.2.fq.gz \
    > ${uuid}.aligned.sam
    """
}

//remove duplicates using samtools v 1.9
process removeDuplicates{

    input:
    	set uuid, file("${uuid}.aligned.sam") from bwa_mapped
    	file "*" from ref_index

    output:
    	set uuid, file("${uuid}.bam"), file("${uuid}.bam.bai") into dup_removed
    
    cpus threads
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


//run samtools mpileup - creates BCF containing genotype likelihoods 
process mpileup{

    input:
    	set uuid, file("${uuid}.bam"), file("${uuid}.bam.bai") from dup_removed
    	file refFasta
    	file "*" from ref_index
 
    output:
    	set uuid, file("pileup.bcf") into pileup
   
   	cpus 1
    tag "${getShortId(uuid)}"

	//use bcftools mpileup to generate vcf file
	//mpileup genearates the likelihood of each base at each site
 	"""
    bcftools mpileup -Ou -f $refFasta ${uuid}.bam > pileup.bcf
    """

}


//call SNPs using samtools call from mpileup file
process snpCall{

    input:
    	set uuid, file("pileup.bcf") from pileup
    	file refFasta
    	file "*" from ref_index
 
    output:
    	set uuid, file("${uuid}.bcf"), file("${uuid}.allsites.bcf") into snps_called
   
    cpus 1
    tag "${getShortId(uuid)}"

	//call converts pileup to actual variants in the BCF or VCF file
	//norm, normalises indels
		//-m +any option to join biallelic sites into multiallelic records
		// and do this for any (i.e. SNPs and indels)
	
    """   
    # call variants only
    # 	-m use multiallelic model
    # 	-v output variants only
    bcftools call -Ou -m -v --ploidy 1 pileup.bcf | \
    	bcftools norm -f $refFasta -m +any -Ou -o ${uuid}.bcf
    	
    # call all sites
    bcftools call -Ou -m --ploidy 1 pileup.bcf | \
    	bcftools norm -f $refFasta -m +any -Ou -o ${uuid}.allsites.bcf
    """

}


//filter SNPS
process filterSnps{

    input:
    	set uuid, file("${uuid}.bcf"), file("${uuid}.allsites.bcf") from snps_called
    	file refFasta
    	file "*" from ref_index
 
    output:
    	set uuid, file("${uuid}.snps.vcf.gz"), file("${uuid}.snps.vcf.gz.csi"), 
    		file("${uuid}.zero_coverage.vcf.gz"), file("${uuid}.zero_coverage.vcf.gz.csi") into filtered_snps
   
    cpus 1
    tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.snps.*"
	publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.indels.*"
	publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.zero_coverage.*"

	//use bcftools to filter normalised bcf file of variants from pileup and call
	//use one line for each filter condition and label
	//create index at end for random access and consensus calling
	
	//filters
		// quality >30
		// one read in each direction to support variant
		// not in a repeat region
		// consensus of >75% reads to support alternative allele
		// mask SNPs within 3 bp of INDEL
		// require high quality depth of 5 for call
	
    """
    #annotate vcf file with repetitive regions
	bcftools annotate -a ${refFasta.baseName}.rpt_mask.gz -c CHROM,FROM,TO,RPT \
		-h ${refFasta.baseName}.rpt_mask.hdr ${uuid}.bcf -Ou -o ${uuid}.masked.bcf
    
    #filter vcf
    bcftools filter -S . -s Q30 -e '%QUAL<30' -Ou ${uuid}.masked.bcf | \
    	bcftools filter -S . -s OneEachWay -e 'DP4[2] == 0 || DP4[3] ==0' -m+ -Ou | \
    	bcftools filter -S . -s RptRegion -e 'RPT=1' -m+ -Ou | \
    	bcftools filter -S . -s Consensus75 -e '((DP4[2]+DP4[3])/(DP4[0]+DP4[1]))<3' -m+ -Ou | \
    	bcftools filter -S . -s SnpIndel --SnpGap 3 -m+ -Ou | \
    	bcftools filter -S . -s HQDepth5 -e '(DP4[2]+DP4[3])<=5' -m+ -Oz -o ${uuid}.all.vcf.gz
    
    #create vcf file with just SNPs
    bcftools filter -i 'TYPE="snp"' -m+ -Oz -o ${uuid}.snps.vcf.gz ${uuid}.all.vcf.gz 
    bcftools index ${uuid}.snps.vcf.gz
    
    #create vcf file with just INDELs
    bcftools filter -i 'TYPE!="snp"' -m+ -Oz -o ${uuid}.indels.vcf.gz ${uuid}.all.vcf.gz 
    bcftools index ${uuid}.indels.vcf.gz
    
    #create vcf file with just zero depth sites
    bcftools filter -e 'DP>0' -Oz -o ${uuid}.zero_coverage.vcf.gz ${uuid}.allsites.bcf
    bcftools index ${uuid}.zero_coverage.vcf.gz
    """
    

}

//generate consensus fasta file
process consensusFa{

	input:
		set uuid, file("${uuid}.snps.vcf.gz"), file("${uuid}.snps.vcf.gz.csi"), 
    		file("${uuid}.zero_coverage.vcf.gz"), file("${uuid}.zero_coverage.vcf.gz.csi") from filtered_snps
		file refFasta
	
	output:
		set uuid, file("${uuid}.fa") into fa_file
	
	cpus 1
	tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/fasta", mode: 'copy', pattern: "${uuid}.*"

	// call consensus sequence
		// -S flag in bcftools filter sets GT (genotype) to missing, with -M flag here
		// setting value to N
	"""
	cat $refFasta | bcftools consensus -H 1 -M "N" ${uuid}.snps.vcf.gz > tmp.fa
	samtools faidx tmp.fa 
	cat tmp.fa | bcftools consensus -H 1 -M "-" ${uuid}.zero_coverage.vcf.gz > ${uuid}.fa
	"""
}


