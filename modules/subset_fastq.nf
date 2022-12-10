process SUBSET_FASTQ {
    tag "Subsample FASTQ files for ${sample_id}"
    publishDir "${params.outdir}/subset_fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)
    val n_reads

    output:
    tuple val(sample_id), path("out/*fastq.gz"), emit: fq 
    path "out/logs/slurm*log", emit: log
    path "out/logs/version.txt", emit: version

    script:
    """
    fqsub.sh \
        --R1_in ${fq_pair[0]} \
        --n_reads ${n_reads} \
        --outdir out
    
    cp .command.log out/logs/slurm-${sample_id}.log
    """
}
