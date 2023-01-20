process MERGE_BAM {
    publishDir "${params.outdir}/map2genome", mode: 'copy'

    input:
    path bamfiles

    output:
    path "merged.bam", emit: bam
    path "logs/slurm*log", emit: log
    path "logs/*_version.txt", emit: version
    
    script:
    """
    merge_bam.sh \
        --outfile merged.bam \
        ${bamfiles}
    
    cp .command.log logs/slurm-merge.log
    """
}
