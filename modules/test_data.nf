process get_fasta {
  	label 'python3'
	output:

	tuple val("test"), path("test.fasta")

	script:
	"""
	wget --no-check-certificate -q  -O - "https://osf.io/download/cv543/" > test.fasta
	"""
    stub:
    """
    touch test.fasta
    """
}

process get_fastq{
  	label 'python3'
	output:
	tuple val("test"), path("test.fastq.gz")
	script:
	"""
	wget --no-check-certificate -q  -O - "https://osf.io/download/68b828e61daee02e98bf5c1b/" > test.fastq.gz
	"""
    stub:
    """
    touch test.fastq.gz
    """    
}

