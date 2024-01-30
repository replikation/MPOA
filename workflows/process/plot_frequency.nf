process plot_frequency {
        label 'ggplot2'
        publishDir "${params.output}", mode: 'copy'
        errorStrategy 'ignore'
    input:
        path(frequency)
    output:
        path("frequency.svg"), optional: true

    script:
        """
        cat ${frequency} | grep -v "nan" > input.tsv
        violin_chart.R
        """

        stub:
        """
        touch frequency.svg
        """
}


