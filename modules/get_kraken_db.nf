process GET_KRAKEN_DB {
    tag "Download a Kraken database"
    publishDir "${params.outdir}/kraken", mode: 'copy'

    input:
    val kraken_db_url

    output:
    path "kraken_db"
    
    script:
    """
    mkdir -p kraken_db
    cd kraken_db

    wget -q -O kraken_db.tar.gz ${kraken_db_url}
    tar -xzvf kraken_db.tar.gz
    """
}
