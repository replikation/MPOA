#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/*
* Nextflow -- mask regions for outbreak analysis
* Author: christian.jena@gmail.com
*/


/************************** 
* Help messages & user inputs & checks
**************************/
// help message
    if (params.help) { exit 0, helpMSG() }


// error codes
    if (params.profile) {
        exit 1, "--profile is WRONG use -profile" }

    if ( workflow.profile == 'standard' ) { 
        "NO EXECUTION PROFILE SELECTED, using [-profile local,docker]" }


    if ( !params.fastq &&  !params.fasta ) {
        exit 1, "input missing, use [--fasta] and [--fastq]"}

 
// nanopore reads
    if (params.fastq) { nano_input_ch = Channel
            .fromPath( params.fastq, checkIfExists: true)
            .map { file -> tuple(file.simpleName, file) }
            }

    if (params.fasta) { fasta_input_ch = Channel
            .fromPath( params.fasta, checkIfExists: true)
            .map { file -> tuple(file.simpleName, file) }
            }





/************************** 
* WORKFLOWS
**************************/

include { mask_regions_wf } from './workflows/mask_regions.nf' 

/************************** 
* MAIN WORKFLOW 
**************************/

workflow {
defaultMSG()

mask_regions_wf(nano_input_ch, fasta_input_ch)
                 
}



/*************  
* LOG INFO DEFINITIONS
*************/

def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________
    
    Workflow: MPOA - Mask Pathogens for Outbreak Analysis

    ${c_yellow}Usage example:${c_reset}
    nextflow run nanozoo/MPOA --fastq '*.fastq' --fasta '*.fasta'

    ${c_yellow}Inputs (Mandatory):
     ${c_green}--fastq ${c_reset}        e.g.: 'sample1.fastq' or '*.fastq' or '*/*.fastq'
     ${c_green}--fasta ${c_reset}        e.g.: 'sample1.fasta' or '*.fasta' 
    Reads and Genome files are mateched based on the first word before the first dot in their filename.
      Matches: ${c_green}Sample1${c_reset}.clean.fasta ${c_green}Sample1${c_reset}.fastq.gz
      This NOT: ${c_yellow}clean${c_reset}.Sample1.clean.fasta ${c_yellow}Sample1${c_reset}.fastq.gz

${c_yellow}Options  (optional)${c_reset}
     --cores         amount of cores for a process (local use) [default: $params.cores]
     --max_cores     max amount of cores for poreCov to use (local use) [default: $params.max_cores]
     --memory        available memory [default: $params.memory]
     --output        name of the result folder [default: $params.output]
     --workdir       defines the path where nextflow writes tmp files [default: $params.workdir]

    ${c_yellow}Profiles:${c_reset}
     -profile               local,docker
                             ${c_reset}

    """.stripIndent()
}

def defaultMSG() {
    log.info """
    MPOA
    \u001B[1;30m______________________________________\033[0m

    \u001B[36mMask Pathogens for Outbreak Analysis\033[0m
    \u001B[1;30m______________________________________\033[0m
    Profile:                $workflow.profile
    Current User:           $workflow.userName
    Nextflow-version:       $nextflow.version
    Starting time:          $nextflow.timestamp
    Workflow hash:          $workflow.commitId
    Workflow revision:      $workflow.revision

        --workdir           $params.workdir
        --output            $params.output
        --cores             $params.cores
        --max_cores         $params.max_cores
        --mem               $params.mem
    \u001B[1;30m______________________________________\033[0m
    """.stripIndent()
}

