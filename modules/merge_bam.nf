process MERGE_BAM {
    tag "Merge BAM files for genome-guided assembly"
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
        -o merged.bam \
        ${bamfiles}
    
    cp .command.log logs/slurm-merge.log
    """
}
