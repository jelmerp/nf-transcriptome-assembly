process RCORRECTOR {
    tag "Rcorrector read correction"
    publishDir "${params.outdir}/rcorrector", mode: 'copy'

    input:
    val fq_list

    output:
    path "*fq.gz", emit: fq
    path "logs/slurm*log", emit: log
    path "fofn*txt", emit: fofn

    script:
    """
    echo "$fq_list" | grep . > fofn.txt
    subset_id=\$(head -n 1 fofn.txt | xargs -I{} basename {} .fastq.gz)

    rcorrector.sh \
        --fofn fofn.txt \
        --outdir .
    
    mv fofn.txt fofn_\${subset_id}.txt
    cp .command.log logs/slurm-rcorrector_\${subset_id}.log
    """
}
