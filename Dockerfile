################## BASE IMAGE ###################### 
FROM continuumio/miniconda3:latest  

################## METADATA ###################### 
LABEL about.summary="Test Dockerfile for running pipelines"

################## MAINTAINER ###################### 
MAINTAINER David Eyre <david.eyre@bdi.ox.ac.uk>  

################## INSTALLATION ######################  
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install bwa=0.7.15
RUN conda install fastqc=0.11.8
WORKDIR /data/