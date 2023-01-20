process MAP2TRANSCRIPTOME {
    tag "Map the reads back to the transcriptome"

    input:
    path assembly
    path fq_fofn

    output:
    path "map2trans.bam", emit: bam
    path "logs/slurm.log", emit: log
    path "logs/version.txt", emit: version
    
    script:
    """
    bowtie2.sh \
        --assembly ${assembly} \
        --fq_fofn ${fq_fofn} \
        --bam map2trans.bam

    cp .command.log logs/slurm.log
    """
}
