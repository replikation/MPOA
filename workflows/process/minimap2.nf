process minimap2 {
        label 'minimap2'
        publishDir "${params.output}/${name}/", mode: 'copy'
    input:
        tuple val(name), path(fasta), path(reads)
  	output:
    	tuple val(name), file("${name}.masked.fasta"), env(MASKREGIONS)
  	script:
    """
    minimap2 -t ${task.cpus} -o ${name}.sam -ax map-ont ${fasta} ${reads}
    samtools view -bS ${name}.sam | samtools sort - -@ ${task.cpus} -o ${name}.minimap.sorted.bam

    # consensus
    samtools consensus -@ ${task.cpus} -f fasta -X r10.4_sup ${name}.minimap.sorted.bam -o ${name}.masked.fasta

    # get Masked Bases
    MASKREGIONS=\$(grep -v ">" ${name}.masked.fasta | grep -o "N" | wc -l)

    """
    stub:
    """
    touch ${name}_masked.fasta
    MASKREGIONS="500"
    """
}