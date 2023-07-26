process get_frequency {
        label 'samtools'
        publishDir "${params.output}/${name}/pileup/", mode: 'copy'
    input:
        tuple val(name), path(fasta), path(positions), path(bam), path(reference)

  	output:
    	tuple val(name), file("${name}_Frequency.tsv"), path("${name}_pileup.tsv")
  	script:
    """
    # create subbam base on fasta name(?)

    # remove headers and create one correct bed position file for samtools (name \t position)
    tail -q -n+2 ${positions}  | awk -v OFS='\\t' '{print \$1,\$5}' > clean.bed

    # pileup for all positions
    samtools index ${bam}
    samtools mpileup ${bam} -f ${fasta} --positions clean.bed --output ${name}_pileup.tsv

    # count stuff per degenerate base
    ## Y
    YT_COUNT=\$(awk -F'\\t' '{if(\$3=="Y")print \$5}' ${name}_pileup.tsv | sed 's/[^tT]//g' | awk '{ print length }')
    Y_TOTAL=\$(awk -F'\\t' '{if(\$3=="Y")print \$5}' ${name}_pileup.tsv | sed 's/[^cCtT]//g' | awk '{ print length }')

    paste -d'\\t' <(echo "\$Y_TOTAL") <(echo "\$YT_COUNT") | awk -v OFS='\\t' '{print \$2 / \$1}'| sed 's/^/T\\t/' | sed 's/^/Y\\t/' >> ${name}_Frequency.tsv

    ## R
    RA_COUNT=\$(awk -F'\\t' '{if(\$3=="R")print \$5}' ${name}_pileup.tsv | sed 's/[^aA]//g' | awk '{ print length }')
    R_TOTAL=\$(awk -F'\\t' '{if(\$3=="R")print \$5}' ${name}_pileup.tsv | sed 's/[^aAgG]//g' | awk '{ print length }')
    
    paste -d'\\t' <(echo "\$R_TOTAL") <(echo "\$RA_COUNT") | awk -v OFS='\\t' '{print \$2 / \$1}'| sed 's/^/A\\t/' | sed 's/^/R\\t/' >> ${name}_Frequency.tsv

    ## W
    WA_COUNT=\$(awk -F'\\t' '{if(\$3=="W")print \$5}' ${name}_pileup.tsv | sed 's/[^aA]//g' | awk '{ print length }')
    W_TOTAL=\$(awk -F'\\t' '{if(\$3=="W")print \$5}' ${name}_pileup.tsv | sed 's/[^aAtT]//g' | awk '{ print length }')
    
    paste -d'\\t' <(echo "\$W_TOTAL") <(echo "\$WA_COUNT") | awk -v OFS='\\t' '{print \$2 / \$1}'| sed 's/^/A\\t/' | sed 's/^/W\\t/' >> ${name}_Frequency.tsv

    ## S
    SG_COUNT=\$(awk -F'\\t' '{if(\$3=="S")print \$5}' ${name}_pileup.tsv | sed 's/[^gG]//g' | awk '{ print length }')
    S_TOTAL=\$(awk -F'\\t' '{if(\$3=="S")print \$5}' ${name}_pileup.tsv | sed 's/[^cCgG]//g' | awk '{ print length }')
    
    paste -d'\\t' <(echo "\$S_TOTAL") <(echo "\$SG_COUNT") | awk -v OFS='\\t' '{print \$2 / \$1}'| sed 's/^/G\\t/' | sed 's/^/S\\t/' >> ${name}_Frequency.tsv

    ## M
    MA_COUNT=\$(awk -F'\\t' '{if(\$3=="M")print \$5}' ${name}_pileup.tsv | sed 's/[^aA]//g' | awk '{ print length }')
    M_TOTAL=\$(awk -F'\\t' '{if(\$3=="M")print \$5}' ${name}_pileup.tsv | sed 's/[^aAcC]//g' | awk '{ print length }')
    
    paste -d'\\t' <(echo "\$M_TOTAL") <(echo "\$MA_COUNT") | awk -v OFS='\t' '{print \$2 / \$1}'| sed 's/^/A\\t/' | sed 's/^/M\\t/' >> ${name}_Frequency.tsv

    ## K
    KG_COUNT=\$(awk -F'\\t' '{if(\$3=="K")print \$5}' ${name}_pileup.tsv | sed 's/[^gG]//g' | awk '{ print length }')
    K_TOTAL=\$(awk -F'\\t' '{if(\$3=="K")print \$5}' ${name}_pileup.tsv | sed 's/[^tTgG]//g' | awk '{ print length }')
    
    paste -d'\\t' <(echo "\$K_TOTAL") <(echo "\$KG_COUNT") | awk -v OFS='\t' '{print \$2 / \$1}'| sed 's/^/G\\t/' | sed 's/^/K\\t/' >> ${name}_Frequency.tsv

    """
}