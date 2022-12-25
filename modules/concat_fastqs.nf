process CONCAT_FASTQS {
    tag "Concatenate all FASTQ files"

    input:
    path fastqs

    output:
    path "*fastq.gz", emit: fq
    path "logs/slurm.log", emit: log
    path "logs/version.txt", emit: version

    script:
    """
    fqcat.sh --outdir . ${fastqs}

    cp .command.log logs/slurm.log
    """
}
