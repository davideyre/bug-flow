################## BASE IMAGE ###################### 
FROM continuumio/miniconda3:latest  

################## METADATA ###################### 
LABEL about.summary="Test Dockerfile for running pipelines"

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

RUN conda install r-base \
				  bwa=0.7.15 \
				  fastqc=0.11.8 \
				  velvet=1.2.10 \
				  spades=3.13.0 \
				  bbmap=38.22 \
				  unicycler=0.4.7 \
				  samtools=1.9 \
				  blast=2.5.0 \
				  perl-velvetoptimiser=v2.2.6

WORKDIR /data/

