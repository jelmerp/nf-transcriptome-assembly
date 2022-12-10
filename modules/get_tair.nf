process GET_TAIR {
    tag "Download the TAIR protein database"
    publishDir "${params.outdir}/dbs/TAIR", mode: 'copy'

    output:
    path "TAIR10_pep_20101214.fasta"

    script:
    """
    wget -O TAIR10_pep_20101214.fasta \
        https://www.arabidopsis.org/download_files/Sequences/TAIR10_blastsets/TAIR10_pep_20101214
    """
}
