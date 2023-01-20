process TRANSABYSS {
    tag "Trans-ABySS - kmer $kmer_size - norm. $norm"
    publishDir "${params.outdir}/transabyss", mode: 'copy'

    input:
    val fq_list
    each kmer_size
    val norm

    output:
    path "transabyss_norm*.fasta", emit: assembly
    path "logs/slurm*log", emit: log
    path "logs/version.txt", emit: version
    
    script:
    """
    echo "$fq_list" | tr "," "\n" | tr -d " ][" | grep . > fofn.txt
    subset_id=\$(head -n 1 fofn.txt | xargs -I{} basename {} .fastq.gz)

    transabyss.sh \
        --fofn fofn.txt \
        --outdir . \
        --id transabyss_norm${norm}_k${kmer_size}_subset\${subset_id} \
        --min_contig_length ${params.min_contig_length} \
        --kmer_size ${kmer_size} \
        --strandedness ${params.strandedness}
    
    mv fofn.txt fofn_norm${norm}_k${kmer_size}_\${subset_id}.txt

    cp .command.log logs/slurm_norm${norm}_k${kmer_size}_\${subset_id}.log
    """
}
