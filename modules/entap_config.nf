process ENTAP_CONFIG {
    tag "Configure EnTap"
    publishDir "${params.outdir}/entap", mode: 'copy'

    input:
    path config
    path protein_dbs

    output:
    path "db_dir/diamond_from_fasta", emit: db_dir
    path "entap_config_final.ini", emit: config
    
    script:
    """
    entap_config.sh \
        --config-in ${config} \
        --config-out entap_config_final.ini \
        --db-dir db_dir \
        --taxon ${params.entap_taxon} \
        --contam ${params.entap_contam} \
        --qcoverage ${params.entap_qcov} \
        --tcoverage ${params.entap_tcov} \
        --evalue ${params.entap_eval} \
        --fpkm ${params.entap_fpkm} \
        ${protein_dbs}

    cp .command.log db_dir/logs/slurm.log
    """
}
