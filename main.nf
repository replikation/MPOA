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
include { checkForWhitespace    } from './lib/functions.nf'
include { helpMSG               } from './lib/messages.nf'
include { defaultMSG            } from './lib/messages.nf'
include { keepMSG               } from './lib/messages.nf'

workflow {
    main:
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
        checkForWhitespace('--fastq', params.fastq)

// nanopore assembly (fasta)
    fasta_input_ch = params.fasta ?  
        Channel.fromPath(params.fasta, checkIfExists: true)
            .map { it -> tuple(it.simpleName, it) } :
        Channel.empty()
        checkForWhitespace('--fasta', params.fasta)

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
        mask_out = mask_regions_wf(fastq_ch, fasta_ch)
        logo_out = identify_logo_wf(mask_out.chromosome_degen, mask_out.bam)

    publish:
        masked_fasta = mask_out.mapper_fasta
        masked_bam = mask_out.bam
        masked_dm = mask_out.masked_dm
        chromosomes = mask_out.chromosome_degen
        plasmids = mask_out.plasmids
        unclassified = mask_out.unclassified
        masked_summary = mask_out.summary
        snp_positions = logo_out.positions.filter { _name, _fasta, positions -> positions }
        snp_regions = logo_out.regions.filter { _name, regions -> regions }
        motif_logos = logo_out.logos
        frequency = logo_out.frequency
        pileup = logo_out.pileup
}

/************************** 
* WORKFLOW OUTPUTS
**************************/
output {
    masked_fasta {path { name  -> "${name[0]}/" }}
    masked_bam {path { name -> "${name[0]}/" }}
    masked_dm {path { name -> "${name[0]}/" }}
    chromosomes {path { name -> "${name[0]}/" }}
    plasmids {path { name-> "${name[0]}/" }}
    unclassified {path { name -> "${name[0]}/" }}
    masked_summary {path '.'}
    snp_positions {path { name -> "${name[0]}/SNP_positions/" }}
    snp_regions {path { name -> "${name[0]}/SNP_regions/" }}
    motif_logos {path 'motif_logos'}
    frequency {path '.'}
    pileup {path { name -> "${name[0]}/pileup/" }}
}