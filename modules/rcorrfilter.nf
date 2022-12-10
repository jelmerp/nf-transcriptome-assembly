process RCORRFILTER {
    tag "Rcorrfilter removal of unfixable reads for ${sample_id}"
    publishDir "${params.outdir}/rcorrfilter", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fq_pair)

    output:
    tuple val(sample_id), path("*fastq.gz"), emit: fq 
    path "logs/slurm*log", emit: log

    script:
    """
    rcorrfilter.sh \
        --R1 ${fq_pair[0]} \
        --outdir .
    
    cp .command.log logs/slurm-${sample_id}.log
    """
}
