process plot_motifs {
        label 'ggplot2'
        publishDir "${params.output}/${name}/", mode: 'copy'
    input:
        //tuple val(name), path(fasta), val(allmasked) 
    output:
        //tuple val(name), path("")

    script:
        """
        # remove motifs with below e.g. 5 occurences
        # write it in the title (n > 4)

        # posit R studio to create it
        """
        stub:
        """
        
        """
}


