process TRIMGALORE {
    tag "TrimGalore to remove adapters and quality-trim for ${sample_id}"
    publishDir "${params.outdir}/trimgalore", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)

    output:
    path "trimmed/*fastq.gz", emit: fq_trimmed
    path "logs/*_trimming_report.txt", emit: trim_report
    path "logs/slurm*log", emit: log
    path "logs/version.txt", emit: version
    path "fastqc/*zip", emit: fqc_zip
    path "fastqc/*html", emit: fqc_html

    script:
    """
    nextseq_arg=""
    [[ ${params.trim_nextseq} = true ]] && nextseq_arg="--nextseq"

    trimgalore.sh \
        --R1 ${fq_pair[0]} \
        --outdir . \
        --quality ${params.trim_qual} \
        --length ${params.trim_len} \
        \${nextseq_arg}
    
    cp .command.log logs/slurm-${sample_id}.log
    """
}
