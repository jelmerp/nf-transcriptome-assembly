process INDEX_GENOME {
    tag "Index a reference genome with STAR"
    publishDir "${params.outdir}/index_genome", mode: 'copy'

    input:
    path ref_fasta

    output:
    path "index_dir", emit: index
    path "index_dir/logs/version.txt", emit: version
    
    script:
    """
    star_index.sh \
        -i ${ref_fasta} \
        -o index_dir
    
    cp .command.log index_dir/logs/slurm.log
    """
}
