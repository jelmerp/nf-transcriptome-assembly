process KRAKEN {
    tag "Kraken contamination detection for ${sample_id}"
    publishDir "${params.outdir}/kraken", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)
    path kraken_db_dir

    output:
    tuple val(sample_id), path("unclassified/*fastq.gz"), emit: fq //No longer used
    path "unclassified/*fastq.gz", emit: fq_list
    path "classified/*fastq.gz", emit: fq_classified
    tuple val(sample_id), path("*main.txt"), emit: main
    path "*report.txt", emit: report
    path "logs/slurm*log", emit: log
    path "logs/version.txt", emit: version
    
    script:
    """
    kraken.sh \
        --infile ${fq_pair[0]} \
        --outdir . \
        --db-dir ${kraken_db_dir} \
        --classified-out \
        --unclassified-out \
        --confidence 0.5
    
    cp .command.log logs/slurm-${sample_id}.log
    """
}
