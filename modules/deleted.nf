
//include { DOWNLOAD_SWISSPROT; DOWNLOAD_EGGNOG_SQL; DOWNLOAD_EGGNOG_DIAMOND }

// For config:
//     withName: 'DOWNLOAD_EGGNOG_DIAMOND' {
//        conda = "/fs/ess/PAS0471/jelmer/conda/diamond"
//    }

WORKFLOW {
    if (swissprot_db == false) {
        swissprot_ch = DOWNLOAD_SWISSPROT()
    } else {
        swissprot_ch = Channel.fromPath(swissprot_db) 
    }
    ref_dbs_ch = refseq_ch.mix(nr_ch, swissprot_ch).collect()
    
    eggnog_sql_ch = DOWNLOAD_EGGNOG_SQL()
    eggnog_diamond_ch = DOWNLOAD_EGGNOG_DIAMOND()
}

process DOWNLOAD_SWISSPROT {
    tag "Download the SwissProt database"
    publishDir "${params.outdir}/dbs/swissprot", mode: 'copy'

    output:
    path "uniprot_sprot.fasta"

    script:
    """
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    gunzip uniprot_sprot.fasta.gz
    """
}

process DOWNLOAD_EGGNOG_SQL {
    tag "Download the EggNOG SQL database"
    publishDir "${params.outdir}/dbs/eggnog", mode: 'copy'

    output:
    path "eggnog.db"

    script:
    """
    wget http://eggnog5.embl.de/download/eggnog_4.1/eggnog-mapper-data/eggnog.db.gz
    gunzip eggnog.db.gz
    """
}

process DOWNLOAD_EGGNOG_DIAMOND {
    tag "Download the EggNOG DIAMOND database"
    publishDir "${params.outdir}/dbs/eggnog", mode: 'copy'

    output:
    path "eggnog_proteins.dmnd"

    script:
    """
    wget http://eggnog5.embl.de/download/eggnog_4.1/eggnog-mapper-data/eggnog4.clustered_proteins.fa.gz

    diamond makedb \
        --threads ${task.cpus} \
        --in eggnog4.clustered_proteins.fa.gz \
        --db eggnog_proteins
    """
}
