process ENTAP_CONFIG {
    tag "Configure EnTap"
    publishDir "${params.outdir}/entap", mode: 'copy'

    input:
    path config

    output:
    path "entap_config_final.ini", emit: config
    path "logs/slurm.log", emit: log

    script:
    """
    entap_config.sh \
        --config_in ${config} \
        --config_out "entap_config_final.ini" \
        --db_dir db_dir \
        --taxon ${params.entap_taxon} \
        --contam ${params.entap_contam} \
        --qcoverage ${params.entap_qcov} \
        --tcoverage ${params.entap_tcov} \
        --evalue ${params.entap_eval} \
        --fpkm ${params.entap_fpkm}

    cp .command.log logs/slurm.log
    """
}
