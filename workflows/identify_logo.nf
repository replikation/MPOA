include { count_positions } from './process/count_positions.nf'
include { extract_sequence } from './process/extract_sequence.nf'
include { plot_motifs } from './process/plot_motifs.nf'
include { get_frequency } from './process/get_frequency.nf'
include { plot_frequency} from './process/plot_frequency.nf'

workflow identify_logo_wf {
    take:   
        fasta // masked and depth masked
        bam // masked bam (not dm masked)
    main: 

        // get positions
        positions_ch = count_positions(fasta)

        // extract regions
        regions_ch = extract_sequence(positions_ch)

        // visualize motifs
        logos_ch = plot_motifs(regions_ch)

        frequency_ch = channel.empty()
        pileup_ch = channel.empty()

        // plot frequency per deg. base
        if (params.frequency) {
            pileup_ch = get_frequency(
                positions_ch // val path( fastadm and masked), path (positions deg bases)
                .join(bam)
            )

            // collect all frequency results
            frequency_ch = plot_frequency(pileup_ch.map{it -> it[1]}.collect())
        }      

    emit:
        positions = positions_ch
        regions = regions_ch
        logos = logos_ch
        frequency = frequency_ch
        pileup = pileup_ch
}