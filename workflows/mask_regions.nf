include { minimap2 } from './process/minimap2.nf'

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

    // report masked status
    report_ch = minimap2.out.view { name, file, masked -> "Sample: $name has $masked masked bases" }

    emit: minimap2.out


}