process plasflow {
        label 'plasflow'
        publishDir "${params.output}/${name}/", mode: 'copy'
    input:
        tuple val(name), path(fasta), val(allmasked) 
    output:
        tuple val(name), path("${name}_chromosomes.fasta"), emit: chromosomes  optional true
        tuple val(name), path("${name}_plasmids.fasta"), emit: plasmids  optional true
        tuple val(name), path("${name}_unclassified.fasta"), emit: unclassified  optional true
        tuple val(name), val(allmasked), env(MASKREGIONS), emit: report
    script:
        """
        PlasFlow.py --input ${fasta} --output ${name} --threshold 0.7

        # remove empty file
        find . -name "*.fasta" -type f -size 0 -print0 | xargs -0 echo rm

        # get Masked Bases
        MASKREGIONS=\$(grep -v ">" ${name}_chromosomes.fasta | grep -o "N" | wc -l)
        """
        stub:
        """
        touch  ${name}_chromosomes.fasta ${name}_plasmids.fasta ${name}_unclassified.fasta
        """
}