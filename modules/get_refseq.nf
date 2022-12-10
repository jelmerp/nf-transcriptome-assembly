process GET_REFSEQ {
    tag "Download an NCBI RefSeq database"
    publishDir "${params.outdir}/dbs/refseq", mode: 'copy'

    input:
    val refseq_type

    output:
    path "refseq_*.fasta"
    
    script:
    """
    wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/${refseq_type}/*protein*faa.gz
    cat complete* | gunzip -c > refseq_${refseq_type}.fasta
    """
}
