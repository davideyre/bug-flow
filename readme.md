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
docker build -t bug-flow .
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
From within the locally cloned repository
```
nextflow run bug-flow.nf
```

Alternatively to resume a partially completed run where the intermediate files have been saved:
```
nextflow run bug-flow.nf -resume
```

Periodically you may wish to delete the working files created by the pipeline, the contents of the work/ directory can be safely removed, but any cached parts of the pipeline will be lost

## Notes
If you run the docker container directly and access it interactively this may cause strange behaviour with the pipeline if you use the same tag.  

David Eyre
david.eyre@bdi.ox.ac.uk
23 January 2019