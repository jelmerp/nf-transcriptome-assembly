process TRINITY_GUIDED {
    publishDir "${params.outdir}/trinity_guided", mode: 'copy'

    input:
    path bam

    output:
    path "trinity_out/trinity*fasta", emit: assembly
    path "trinity_out/trinity*gene2trans", emit: gene2trans
    path "trinity_out/logs/slurm*log", emit: log
    path "trinity_out/logs/version.txt", emit: version
    
    script:
    """
    subset_id=\$(basename $bam .bam)

    trinity.sh \
        --bam ${bam} \
        --outdir trinity_out \
        --genome_guided_max_intron 250000 \
        --strandedness ${params.strandedness}
    
    cp -v .command.log trinity_out/logs/slurm_\${subset_id}.log
    """
}
