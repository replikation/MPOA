include { minimap2 } from './process/minimap2.nf'
include { plasflow } from './process/plasflow.nf'

workflow mask_regions_wf {
    take:   
        fastq
        fasta
    main: 

    // join channels
    combined_ch = fasta.join(fastq, by:0)

    // check if everything is alright
    no_match = combined_ch.ifEmpty{ log.info "\033[0;33mCould not match any reads to genomes, please read the help via --help\033[0m" }

    minimap2(combined_ch)
    plasflow(minimap2.out.fasta)

    // report masked status
    report_ch = plasflow.out.report.view { name, all, chr -> "$name: Total masked Bases $all with $chr on the chromosome only." }

    report_write_ch = plasflow.out.report
        .collectFile(seed: 'name,total masked bases,masked bases chromosome\n', 
                    storeDir: params.output + "/") {
                    row -> [ "masked_bases_summary.csv", row[0] + ',' + row[1] + ',' + row[2] + '\n']
                    }

    emit: minimap2.out.fasta


}