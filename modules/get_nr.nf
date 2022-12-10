process GET_NR {
    tag "Download the NCBI NR database"
    publishDir "${params.outdir}/dbs/nr", mode: 'copy'

    output:
    path "nr_database.fasta"

    script:
    """
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
    
    gunzip -c nr.gz > nr_database.fasta
    """
}
