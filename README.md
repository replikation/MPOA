<p align="center">
  <img src="data/logo/mobile_logo.png" width="800" title="Workflow">
</p>

**MPOA | Workflow to mask pathogens for outbreak analysis**   
===
![](https://img.shields.io/github/v/release/replikation/MPOA)
![](https://img.shields.io/badge/uses-Docker-blue.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)


[![Twitter Follow](https://img.shields.io/twitter/follow/maralohde.svg?style=social)](https://twitter.com/maralohde) 
[![Twitter Follow](https://img.shields.io/twitter/follow/gcloudChris.svg?style=social)](https://twitter.com/gcloudChris) 

>https://www.biorxiv.org/content/10.1101/2023.09.15.556300v1

## What is this Repo?
* A quick workflow to mask fasta files for outbreak analysis
* It mitigates false base calls from e.g., modified bases by masking all such positions that are "uncertain"
* Figures and overviews are provided to assess weather your samples might be affected by such an issue

<p align="center">
  <img src="data/figures/figure_4_flowchart.drawio.png" width="500" title="Workflow">
</p>

# Quick installation
## 1.1 Nextflow (the workflow manager)
* MPOA needs [Nextflow](https://www.nextflow.io/index.html) and java run time (default-jre)
    * install java run time via:  `sudo apt install -y default-jre`
    * install Nextflow e.g.::  `curl -s https://get.nextflow.io | bash && sudo mv nextflow /bin && sudo chmod 770 /bin/nextflow`
## 1.2 Docker
* Installation [here](https://docs.docker.com/v17.09/engine/installation/linux/docker-ce/ubuntu/#install-docker-ce) (recommended), alternatively via: `sudo apt install -y docker`
* Add Docker to the user: `sudo usermod -a -G docker $USER`

# Quick start

* Reads and genomes are matched by their name before the first "."
* E.g. genome file **sample1**.test.fasta matches with read file **sample1**.fastq.gz

```bash
# help
nextflow run replikation/MPOA -profile local,docker --help

#example run
nextflow run nanozoo/MPOA --fastq '*.fastq' --fasta '*.fasta' -profile local,docker
```
