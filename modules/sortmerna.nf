process SORTMERNA {
    tag "SortMeRNA rRNA removal for ${sample_id}"
    publishDir "${params.outdir}/sortmerna", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)

    output:
    tuple val(sample_id), path("unmapped/*fastq.gz"), emit: fq_unmapped
    path "mapped/*fastq.gz", emit: fq_mapped 
    path "logs/*log", emit: log
    path "logs/version.txt", emit: version
    
    script:
    """
    sortmerna.sh \
        --R1 ${fq_pair[0]} \
        --outdir .
    
    cp .command.log logs/slurm-${sample_id}.log
    """
}
