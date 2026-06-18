process count_positions {
        label 'seqkit'

    input:
        tuple val(name), path(fasta)
    output:
        tuple val(name), path(fasta), path("*_positions.txt"), optional: true
    script:
        """
        grep -o . <<< "WSMKRYBDHV" | while read BASE; do 
                seqkit locate --ignore-case --only-positive-strand  -r --pattern "\${BASE}" ${fasta} > \${BASE}_positions.txt

                if [ "\$( wc -l <\${BASE}_positions.txt )" -eq 1 ]; then
                    rm -f \${BASE}_positions.txt
                fi
        done
        """
    stub:
        """
        touch *_positions.txt
        """
}

