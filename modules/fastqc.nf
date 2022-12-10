process FASTQC {
    tag "FastQC to QC FASTQ files for ${sample_id}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)

    output:
    path "*zip", emit: zip
    path "*html", emit: html
    path "logs/slurm*log", emit: log

    script:
    """
    fastqc ${fq_pair}
    
    mkdir -p logs
    cp .command.log logs/slurm-${sample_id}.log
    """
}
