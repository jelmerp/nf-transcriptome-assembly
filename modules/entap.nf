process ENTAP {
    tag "Annotate an assembly with EnTap"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path assembly
    path db_dir
    path config
    path bam

    output:
    path "entap", emit: out
    path "logs/slurm*log", emit: log
    
    script:
    """
    entap.sh \
        --assembly ${assembly} \
        --config ${config} \
        --db-dir ${db_dir} \
        --bam ${bam} \
        --outdir entap

    cp .command.log entap/logs/slurm.log
    """
}
