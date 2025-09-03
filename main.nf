#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
* Nextflow -- mask regions for outbreak analysis
* Author: christian.jena@gmail.com
*/

/************************** 
* WORKFLOWS
**************************/
include { get_fasta; get_fastq  } from './modules/test_data.nf'
include { mask_regions_wf       } from './workflows/mask_regions.nf' 
include { identify_logo_wf      } from './workflows/identify_logo.nf'

/* Custom functions and Messages */
include { exit_with_error       } from './lib/functions.nf'
//include { checkForWhitespace    } from './lib/functions.nf'
include { helpMSG               } from './lib/messages.nf'
include { defaultMSG            } from './lib/messages.nf'
include { keepMSG               } from './lib/messages.nf'

workflow {
/**************************
* Help messages & checks
**************************/

// help message
    params.help             ? { exit 0, helpMSG() }()       : defaultMSG()

// error codes
    if (!workflow.profile.contains('test') && ![params.fasta, params.fastq].any{it}) {
                exit_with_error(6, "Input missing See also --help for more.")}
    if (params.profile && !params.profile.startsWith('-'))                  {exit_with_error(1, "[--profile] is WRONG use [-profile].")}
    if (!workflow.profile.contains('test') && ![params.fasta,params.fastq]) {exit_with_error(1, "Provide both --fasta and --fastq input files.")}
    if ( workflow.profile == 'standard' )                                   {exit_with_error(1, "NO EXECUTION PROFILE SELECTED, using [-profile local,docker]")}

    workflow.onError = {
        println workflow.errorReport
    }

 /************************** 
* INPUTs
**************************/

// nanopore reads (fastq)
    fastq_input_ch = params.fastq ?
        Channel.fromPath(params.fastq, checkIfExists: true)
            .map { it -> tuple(it.simpleName, it) } :
        Channel.empty()
        //checkForWhitespace('--fastq', params.fastq)

// nanopore assembly (fasta)
    fasta_input_ch = params.fasta ?  
        Channel.fromPath(params.fasta, checkIfExists: true)
            .map { it -> tuple(it.simpleName, it) } :
        Channel.empty()
        //checkForWhitespace('--fasta', params.fasta)

/************************** 
* WORKFLOW LOGIC
**************************/
    // 0. Test data
    fasta_test_ch = workflow.profile.contains('test') ? get_fasta() : Channel.empty()
    fastq_test_ch = workflow.profile.contains('test') ? get_fastq() : Channel.empty()
    
    // 1. Main Workflow

        fastq_ch = fastq_input_ch.mix(fastq_test_ch)
        fasta_ch = fasta_input_ch.mix(fasta_test_ch)

    if (params.keep_coordinates) {keepMSG()}
        mask_regions_wf(fastq_ch, fasta_ch)
        identify_logo_wf(mask_regions_wf.out.chromosome_degen, mask_regions_wf.out.bam)
}