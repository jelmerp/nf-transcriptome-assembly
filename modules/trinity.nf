process TRINITY {
    tag "Trinity - normalization $norm"
    publishDir "${params.outdir}/trinity", mode: 'copy'

    input:
    val fq_list
    each norm
    // NOTE: The workflow will assume that normalization was done prior to Trinity
    // Therefore, even with 'norm = true', we still pass '--normalize false' to the script

    output:
    path "trinity_out/trinity_norm*fasta", emit: assembly
    path "trinity_out/trinity_norm*gene2trans", emit: gene2trans
    path "trinity_out/fofn_*", emit: fofn
    path "trinity_out/logs/slurm*log", emit: log
    path "trinity_out/logs/version.txt", emit: version
    
    script:
    """
    echo "$fq_list" | tr "," "\n" | tr -d " ][" | grep . > fofn.txt
    subset_id=\$(head -n 1 fofn.txt | xargs -I{} basename {} .fastq.gz)

    trinity.sh \
        --fofn fofn.txt \
        --outfile trinity_out/trinity_norm${norm}_subset\${subset_id}.fasta \
        --normalize false \
        --min_contig_length ${params.min_contig_length} \
        --strandedness ${params.strandedness}
    
    mv -v fofn.txt trinity_out/fofn_norm${norm}_\${subset_id}.txt
    
    cp -v .command.log trinity_out/logs/slurm_norm${norm}_\${subset_id}.log
    """
}
