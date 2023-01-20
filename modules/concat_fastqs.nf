process CONCAT_FASTQS {
    tag "Concatenate all FASTQ files"

    input:
    tuple val(sample_id), path(fq_pair)

    output:
    path "*fastq.gz", emit: fq
    path "logs/slurm.log", emit: log
    path "logs/version.txt", emit: version

    script:
    """
    fqcat.sh \
        --outdir .
        ${fq_pair}

    cp .command.log logs/slurm.log
    """
}
