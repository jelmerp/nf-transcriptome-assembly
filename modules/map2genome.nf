process MAP2GENOME {
    publishDir "${params.outdir}/map2genome", mode: 'copy'

    input:
    val fq_list
    path ref_index

    output:
    path "bam/*bam", emit: bam
    path "fofn*txt", emit: fofn
    path "logs/slurm*log", emit: log
    path "star_logs/*", emit: star_log
    path "logs/version.txt", emit: version
    
    script:
    """
    echo "$fq_list" | tr "," "\n" | tr -d " ][" | grep . > fofn.txt
    subset_id=\$(head -n 1 fofn.txt | xargs -I{} basename {} .fastq.gz)

    star_align.sh \
        --fofn fofn.txt \
        --index_dir ${ref_index} \
        --outdir .

    mv -v fofn.txt fofn_\${subset_id}.txt
    cp -v .command.log logs/slurm-\${subset_id}.log
    """
}
