#!/usr/bin/env nextflow

/* simple example pipeline for mapping and variant calling */

/* TO DO

 - inserts are not saved as have to keep alignment the same length
 - generate rep masks
 
*/

def getShortId( str ) {
  return str.substring(0,8) 
}


// parameters 
params.index = "example_data/file_list.csv"
params.outputPath = "example_output"
params.refFile = "example_data/R00000419.fasta"
params.threads = 6

// initial logging 
log.info "Pipeline Test -- version 0.1"
log.info "Input file              : ${params.index}"
log.info "Output path             : ${params.outputPath}"

// input validation 
outputPath = file(params.outputPath)

Channel
    .fromPath(params.index)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleid, row.uuid, file(row.fq1), file(row.fq2)) }
    .set { samples_ch }

samples_ch.into { samples_ch1; samples_ch2 }

refFasta = file(params.refFile)
threads = params.threads


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


// Map reads to reference genome with BWA MEM
process bwa{

    input:
    	set sampleid, uuid, file(fq1), file(fq2) from samples_ch1
    	file "*" from ref_index
    	file refFasta

    output:
    	set uuid, file("${uuid}.aligned.sam") into bwa_mapped
    
    cpus threads
    tag "${getShortId(uuid)}"

    """
    bwa mem -t ${task.cpus} \
    		$refFasta \
    		$fq1 \
    		$fq2 \
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


