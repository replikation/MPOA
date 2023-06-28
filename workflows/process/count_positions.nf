process count_positions {
        label 'seqkit'
        publishDir "${params.output}/${name}/SNP_positions/", mode: 'copy'
    input:
        tuple val(name), path(fasta)
    output:
        tuple val(name), path(fasta), path("*_positions.txt")
    script:
        """
        grep -o . <<< "WSMKRYBDHV" | while read BASE; do 
                seqkit locate --ignore-case --only-positive-strand  -r --pattern "\${BASE}" ${fasta} > \${BASE}_positions.txt
        done
        """
}

