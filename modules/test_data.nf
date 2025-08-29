process get_fasta {
  	label 'python3'
	output:

	tuple val("test.fasta")

	script:
	"""
	wget --no-check-certificate -q  -O - "https://osf.io/download/p4ny7/" > test.fasta
	"""
    stub:
    """
    touch test.fasta
    """
}

process get_fastq{
  	label 'python3'
	output:
	tuple val("test.fastq.gz")
	script:
	"""
	wget --no-check-certificate -q  -O - "https://osf.io/download/x9m3k/" > test.fastq.gz
	"""
    stub:
    """
    touch test.fastq.gz
    """    
}

