process TRINITY {
    tag "Assembly with Trinity - normalization $norm"
    publishDir "${params.outdir}/trinity", mode: 'copy'

    input:
    path dir_with_all_fqs
    each norm

    output:
    path "trinity_out/trinity_norm*fasta", emit: assembly
    path "trinity_out/trinity_norm*gene_trans_map", emit: gene2trans
    path "trinity_out/logs/slurm*log", emit: log
    path "trinity_out/logs/version.txt", emit: version
    
    script:
    """
    trinity.sh \
        --input ${dir_with_all_fqs} \
        --outdir trinity_out \
        --min_contig_length ${params.min_contig_length} \
        --normalize ${norm} \
        --strandedness ${params.strandedness}
    
    mv -v trinity_out.Trinity.fasta trinity_out/trinity_norm${norm}.fasta
    mv -v trinity_out.Trinity.fasta.gene_trans_map trinity_out/trinity_norm${norm}.gene_trans_map

    cp .command.log trinity_out/logs/slurm.log
    """
}
