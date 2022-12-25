process MAP2GENOME {
    tag "Map reads to a reference genome with STAR for ${sample_id}"
    publishDir "${params.outdir}/map2genome", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)
    path ref_index

    output:
    path "bam/*bam", emit: bam
    path "logs/slurm*log", emit: log
    path "star_logs/*", emit: star_log
    path "logs/version.txt", emit: version
    
    script:
    """
    star_align.sh \
        --R1 ${fq_pair[0]} \
        --index_dir ${ref_index} \
        --outdir .
    
    cp .command.log logs/slurm-${sample_id}.log
    """
}
