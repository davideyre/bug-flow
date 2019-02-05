################## BASE IMAGE ###################### 
FROM continuumio/miniconda3:latest  

################## METADATA ###################### 
LABEL about.summary="BUGflow: nextflow based pipeline for processing bacterial sequencing data"

################## MAINTAINER ###################### 
MAINTAINER David Eyre <david.eyre@bdi.ox.ac.uk>  

################## INSTALLATION ######################  

RUN apt-get clean all && \
 apt-get update && \
 apt-get upgrade -y && \
 apt-get install -y procps && \
 apt-get clean && \
 apt-get purge && \
 rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge


RUN conda install bwa=0.7.15 \
				  fastqc=0.11.8 \
				  spades=3.13.0 \
				  bbmap=38.22 \
				  samtools=1.9 \
				  htslib=1.9 \
				  bcftools=1.9 \
				  blast=2.5.0 \
				  numpy=1.15.4 \
				  biopython=1.72 \
				  openssl=1.0 \
				  gatk4=4.1.0.0 \
				  freebayes=1.2.0 \
				  vcflib=1.0.0_rc2

# need older version of openssl to support samtools 1.9

WORKDIR /data/

