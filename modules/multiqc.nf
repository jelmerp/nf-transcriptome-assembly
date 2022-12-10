process MULTIQC {
    tag "Run MultiQC"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path all_logs
    val file_prefix

    output:
    path "*html", emit: report
    path "logs/slurm*log", emit: log

    script:
    """
    multiqc --filename ${file_prefix}.html --interactive ${all_logs}
    
    mkdir -p logs
    cp .command.log logs/slurm-${file_prefix}.log
    """
}
