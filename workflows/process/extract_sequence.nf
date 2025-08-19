process extract_sequence {
        label 'samtools'
        publishDir "${params.output}/${name}/SNP_regions/", mode: 'copy'
        errorStrategy 'ignore'
    input:
        tuple val(name), path(fasta), path(positions)
    output:
        tuple val(name), path("*.regions.fasta"), optional: true
    script:
        """
        # WSMKRYBDHV
        grep -o . <<< "WSMKRYBDHV" | while read BASE; do 
            awk -v  s=${params.motif} '{print \$1, \$4, \$5-s, \$6+s}' \${BASE}_positions.txt | tail -n+2 | tr "+" ":" | sed 's/ //'| sed 's/ //' | tr " " "-" > \${BASE}_input_samtools.tsv
            
            # if file empty
            if [ -s \${BASE}_input_samtools.tsv ]; then
                    # The file is not-empty.
                    samtools faidx ${fasta} -c -r \${BASE}_input_samtools.tsv > \${BASE}.regions.fasta
            else
                    # The file is empty.
                    echo "skipping samtools for \${BASE}"
            fi

        done
        """
        stub:
        """
        touch ${name}_W.regions.fasta ${name}_S.regions.fasta ${name}_M.regions.fasta ${name}_K.regions.fasta ${name}_R.regions.fasta ${name}_Y.regions.fasta ${name}_B.regions.fasta ${name}_D.regions.fasta ${name}_H.regions.fasta ${name}_V.regions.fasta ${name}_H.regions.fasta
        touch ${name}_masked.fasta ${name}_depth_file.txt ${name}.masked.sorted.bam ${name}.masked.dm.fasta
        BASE_N=5
        BASE_W=7
        BASE_S=1
        BASE_M=3
        BASE_K=0
        BASE_R=0
        """
}