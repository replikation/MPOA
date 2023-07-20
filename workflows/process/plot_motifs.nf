process plot_motifs {
        label 'ggplot2'
        publishDir "${params.output}/motif_logos", mode: 'copy'
        errorStrategy 'ignore'
    input:
        tuple val(name), path(fasta)
    output:
        tuple val(name), path("${name}.svg"), optional: true

    script:
        """
        if [ -f ?.regions.fasta ]; then
    
            create_logo.R ${name}
            mv logo.svg ${name}.svg    
        
        fi

        """

        stub:
        """
        touch ${name}.svg
        """
}


