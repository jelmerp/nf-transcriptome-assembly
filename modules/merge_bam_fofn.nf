process MERGE_BAM {
    publishDir "${params.outdir}/map2genome/merged", mode: 'copy'

    input:
    val bam_list

    output:
    path "*merged.bam", emit: bam
    path "fofn_*", emit: fofn
    path "logs/slurm*log", emit: log
    path "logs/*_version.txt", emit: version
    
    script:
    """
    echo "$bam_list" | tr "," "\n" | tr -d " ][" | grep . > fofn.txt
    subset_id=\$(head -n 1 fofn.txt | xargs -I{} basename {} .bam)

    merge_bam.sh \
        --fofn fofn.txt \
        --outfile \${subset_id}_merged.bam
    
    mv -v fofn.txt fofn_\${subset_id}.txt
    cp .command.log logs/slurm-merge_\${subset_id}.log
    """
}
