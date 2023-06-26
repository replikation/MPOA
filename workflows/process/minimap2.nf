process minimap2 {
        label 'minimap2'
        publishDir "${params.output}/${name}/", mode: 'copy'
    input:
        tuple val(name), path(fasta), path(reads)
  	output:
    	tuple val(name), file("${name}.masked.fasta"), env(MASKREGIONS), emit: fasta
        tuple val(name), file("${name}.masked.sorted.bam"), emit: bam
  	script:
    """
    minimap2 -t ${task.cpus} -o ${name}.sam -ax map-ont ${fasta} ${reads}
    samtools view -bS ${name}.sam | samtools sort - -@ ${task.cpus} -o ${name}.minimap.sorted.bam

    # consensus
    samtools consensus -@ ${task.cpus} --min-depth ${params.depth} -f fasta -X r10.4_sup ${name}.minimap.sorted.bam -o ${name}.masked.fasta
    rm ${name}.minimap.sorted.bam

    # get Masked Bases
    MASKREGIONS=\$(grep -v ">" ${name}.masked.fasta | grep -o "N" | wc -l)

    # rebam to masked file for visuals
    minimap2 -t ${task.cpus} -o ${name}.masked.sam -ax map-ont ${name}.masked.fasta ${reads}
    samtools view -bS ${name}.masked.sam | samtools sort - -@ ${task.cpus} -o ${name}.masked.sorted.bam
    """
    stub:
    """
    touch ${name}_masked.fasta
    MASKREGIONS="500"
    """
}

process minimap2_degen {
        label 'minimap2'
        publishDir "${params.output}/${name}/", mode: 'copy'
    input:
        tuple val(name), path(fasta), path(reads)
  	output:
    	tuple val(name), file("${name}.masked.fasta"), emit: fasta
        tuple val(name), file("${name}.masked.sorted.bam"), emit: bam
        tuple val(name), env(BASE_N), env(BASE_W), env(BASE_S), env(BASE_M), env(BASE_K), env(BASE_R), 
                         env(BASE_Y), env(BASE_B), env(BASE_D), env(BASE_H), env(BASE_V), emit: report
  	script:
    """
    minimap2 -t ${task.cpus} -o ${name}.sam -ax map-ont ${fasta} ${reads}
    samtools view -bS ${name}.sam | samtools sort - -@ ${task.cpus} -o ${name}.minimap.sorted.bam

    # consensus
    samtools consensus -@ ${task.cpus} \
                        --ambig \
                        -f fasta  \
                        --min-depth ${params.depth} \
                        -X r10.4_sup \
                        ${name}.minimap.sorted.bam \
                        -o ${name}.masked.fasta

    # reduce disk footprint                    
    rm ${name}.minimap.sorted.bam

    # get Masked Bases
    BASE_N=\$(grep -v ">" ${name}.masked.fasta | grep -o "N" | wc -l)
    BASE_W=\$(grep -v ">" ${name}.masked.fasta | grep -o "W" | wc -l)
    BASE_S=\$(grep -v ">" ${name}.masked.fasta | grep -o "S" | wc -l)
    BASE_M=\$(grep -v ">" ${name}.masked.fasta | grep -o "M" | wc -l)
    BASE_K=\$(grep -v ">" ${name}.masked.fasta | grep -o "K" | wc -l)
    BASE_R=\$(grep -v ">" ${name}.masked.fasta | grep -o "R" | wc -l)
    BASE_Y=\$(grep -v ">" ${name}.masked.fasta | grep -o "Y" | wc -l)
    BASE_B=\$(grep -v ">" ${name}.masked.fasta | grep -o "B" | wc -l)
    BASE_D=\$(grep -v ">" ${name}.masked.fasta | grep -o "D" | wc -l)
    BASE_H=\$(grep -v ">" ${name}.masked.fasta | grep -o "H" | wc -l)
    BASE_V=\$(grep -v ">" ${name}.masked.fasta | grep -o "V" | wc -l)


    # rebam to masked file for visuals
    minimap2 -t ${task.cpus} -o ${name}.masked.sam -ax map-ont ${name}.masked.fasta ${reads}
    samtools view -bS ${name}.masked.sam | samtools sort - -@ ${task.cpus} -o ${name}.masked.sorted.bam
    """
    stub:
    """
    touch ${name}_masked.fasta ${name}.masked.sorted.bam
    BASE_N=5
    BASE_W=7
    BASE_S=1
    BASE_M=3
    BASE_K=0
    BASE_R=0
    BASE_Y=0
    BASE_B=20
    BASE_D=15
    BASE_H=3
    BASE_V=0
    """
}