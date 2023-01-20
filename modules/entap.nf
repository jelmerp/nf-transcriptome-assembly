process ENTAP {
    tag "Annotate an assembly with EnTap"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path assembly
    path config
    path diamond_dbs
    path bam

    output:
    path "entap", emit: out
    path "logs/slurm*log", emit: log
    
    script:
    """
    entap.sh \
        --assembly ${assembly} \
        --config ${config} \
        --bam ${bam} \
        --outdir entap \
        ${diamond_dbs}

    cp .command.log entap/logs/slurm.log
    """
}
