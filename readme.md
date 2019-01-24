# BUGflow
A bacterial sequencing data pipeline

## Overview
A nextflow-based pipeline for processing bacterial sequencing data generated using Illumina sequencing platforms. Uses a docker image for running the pipeline.

This is intended to be an example of what is possible with nextflow rather than a polished pipeline. However, it does provide a way of approximately replicating the pipeline used by the Modernising Medical Microbiology consortium, based at the University of Oxford

## Installation
Requires a local installation of 
* Docker - https://www.docker.com/get-started
* Java version 8 or later (required for nextflow)
* Nextflow - https://www.nextflow.io

### Clone the repository
Clone the repository locally
```
git clone https://github.com/davideyre/bug-flow.git
```

### Get the docker image
This can be pulled from docker hub
```
docker pull davideyre/bug-flow
```

Alternatively the docker image can be built from the Dockerfile. Within the cloned repository:
```
cd docker
docker build -t davideyre/bug-flow .
```
Note the tag has to match in the `nextflow.config` file.


### Download the example data
From EBI download 2 sets of fastq files
```
cd example_data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR334/008/SRR3349138/SRR3349138_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR334/008/SRR3349138/SRR3349138_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR334/004/SRR3349174/SRR3349174_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR334/004/SRR3349174/SRR3349174_2.fastq.gz
```


## Run the pipeline

### Example pipeline
Having downloaded the example files above to `example_data`, from within the locally cloned repository
```
nextflow run bug-flow.nf
```

Alternatively to resume a partially completed run where the intermediate files have been saved:
```
nextflow run bug-flow.nf -resume
```

### Run your own data
Three input options need to be set
1. `-seqlist my_seqlist.csv`: a CSV file with 4 columns (sampleid, uuid, fq1, fq2). These headers must be used. Sampleid is an arbitrary name for your sample, uuid is a globally unique identifier for your sample that will be used to store output (this can be generated e.g. with `uuidgen` on linux/mac), fq1 and fq2 are the paths to the forward and reverse read fastq files for your samples.
2. `-outputPath my_output_directory`: the path to your output folder
3. `refFile my_ref.fasta`: the path to your reference file to be used for mapping and variant calling.

Examples of each of the three files can be found in `example_data`.

### Nextflow config
By default the workflow will run using up to 6 cores and 16GB of memory, e.g. on a laptop / desktop. An alternative set up can be used. This is specified in the nextflow.config file. More details on potential settings are given here - https://www.nextflow.io/docs/latest/config.html#. Example profiles are provided for running the pipeline on a server (36 cores, 384GB memory) and a SGE cluster. To run the pipeline using the server profile:
```
nextflow run bug-flow.nf -profile server \
	--seqlist example_data/file_list.csv \
	--outputPath example_output \
	--refFile example_data/R00000419.fasta
```

## Notes

### Singularity support
If docker is not available in your setup, it is also possible to create a singularity image. A recipe for building the image is supplied:

```
cd singularity
sudo singularity build bug-flow.img Singularity
```

If you do not have sudo access on the pipeline compute, the image can be built separately and copied to `singularity/bug-flow.img`.

You will the need to modify the nextflow.config file to use Singularity.

Change 

```
process.container = 'davideyre/bug-flow:latest'
docker.enabled = true
```

to

```
process.container = "singularity/bug-flow.img"
singularity.enabled = true
```

Alternatively you can run the pipeline with a profile that uses singularity, e.g. the cluster profile.

### Housekeeping
Periodically you may wish to delete the working files created by the pipeline, the contents of the work/ directory can be safely removed, but any cached parts of the pipeline will be lost

### Notes
If you run the docker container outside of the pipeline and access it interactively this may cause strange behaviour with the pipeline if you use the same image.  

David Eyre
david.eyre@bdi.ox.ac.uk
23 January 2019