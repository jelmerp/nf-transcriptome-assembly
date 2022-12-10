
process DETONATE {
    tag "Evaluate an assembly with Detonate"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path assembly
    path dir_with_all_fqs

    output:
    path "detonate", emit: out
    path "detonate/logs/slurm.log", emit: log
    path "detonate/logs/version.txt", emit: version
    
    script:
    """
    detonate.sh \
        --assembly ${assembly} \
        --fq-dir ${dir_with_all_fqs} \
        --outdir detonate

    cp .command.log detonate/logs/slurm.log
    """
}
