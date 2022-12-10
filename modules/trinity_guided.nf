process TRINITY_GUIDED {
    tag "Assembly with Trinity (genome-guided)"
    publishDir "${params.outdir}/trinity_guided", mode: 'copy'

    input:
    path bam

    output:
    path "trinity_out/trinity_gg.fasta", emit: assembly
    path "trinity_out/trinity_gg.gene_trans_map", emit: gene2trans
    path "trinity_out/logs/slurm*log", emit: log
    path "trinity_out/logs/version.txt", emit: version
    
    script:
    """
    trinity.sh \
        --input ${bam} \
        --outdir trinity_out \
        --genome_guided \
        --genome_guided_max_intron 250000 \
        --strandedness ${params.strandedness}
    
    mv -v trinity_out/Trinity-GG.fasta trinity_out/trinity_gg.fasta
    mv -v trinity_out/Trinity-GG.fasta.gene_trans_map trinity_out/trinity_gg.gene_trans_map

    cp .command.log trinity_out/logs/slurm.log
    """
}
