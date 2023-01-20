process DIAMOND_MAKEDB {
    tag "$fasta"
    publishDir "${params.outdir}/dbs/diamond", mode: 'copy'

    input:
    path fasta

    output:
    path "*dmnd", emit: db
    path "logs/slurm*log", emit: log
    path "logs/*version.txt", emit: version
    
    script:
    """
    diamond_db.sh \
        --infile $fasta \
        --outdir .

    cp .command.log logs/slurm-${fasta}.log
    """
}
