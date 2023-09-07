<p align="center">
  <img src="data/logo/mobile_logo.png" width="800" title="Workflow">
</p>

**MPOA | Workflow to mask pathogens for outbreak analysis**   
===
![](https://img.shields.io/github/v/release/replikation/MPOA)
![](https://img.shields.io/badge/uses-Docker-blue.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)


[![Twitter Follow](https://img.shields.io/twitter/follow/gcloudChris.svg?style=social)](https://twitter.com/gcloudChris) 
[![Twitter Follow](https://img.shields.io/twitter/follow/gcloudChris.svg?style=social)](https://twitter.com/gcloudChris) 



## What is this Repo?
* A quick and simple workflow to mask fasta files for outbreak analysis
* it mitigates false basecalles (e.g. due to modified bases) by masking all such positions that are "uncertain"

# Quick installation
## 1.1 Nextflow (the workflow manager)
* MPOA needs [Nextflow](https://www.nextflow.io/index.html) and java run time (default-jre)
    * install java run time via:  `sudo apt install -y default-jre`
    * install Nextflow e.g.::  `curl -s https://get.nextflow.io | bash && sudo mv nextflow /bin && sudo chmod 770 /bin/nextflow`
## 1.2 Docker
* installation [here](https://docs.docker.com/v17.09/engine/installation/linux/docker-ce/ubuntu/#install-docker-ce) (recommended), alternatively via: `sudo apt install -y docker`
* add Docker to the user: `sudo usermod -a -G docker $USER`

# Quick start

* reads and genomes are matched by their name before the first "."
* e.g. genome file **sample1**.test.fasta matches with read file **sample1**.fastq.gz

```bash
# help
nextflow run replikation/MPOA -profile local,docker --help

#example run
nextflow run nanozoo/MPOA --fastq '*.fastq' --fasta '*.fasta' -profile local,docker
```