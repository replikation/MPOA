include { count_positions } from './process/count_positions.nf'
include { extract_sequence } from './process/extract_sequence.nf'
include { plot_motifs } from './process/plot_motifs.nf'


workflow identify_logo_wf {
    take:   
        fasta
    main: 

        // get positions
        extract_sequence(count_positions(fasta))

        // visualize motifs
        plot_motifs(extract_sequence.out)


    emit:  extract_sequence.out
}