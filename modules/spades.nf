process SPADES {
    tag "SPAdes - kmer $kmer_size - norm. $norm"
    publishDir "${params.outdir}/spades", mode: 'copy'

    input:
    val fq_list
    each kmer_size
    val norm

    output:
    path "spades_norm*.fasta", emit: assembly
    path "logs/slurm*log", emit: log
    path "logs/version.txt", emit: version
    
    script:
    """
    echo "$fq_list" | tr "," "\n" | tr -d " ][" | grep "_R1" > fofn.txt
    subset_id=\$(head -n 1 fofn.txt | xargs -I{} basename {} .fastq.gz)

    spades.sh \
        --fofn fofn.txt \
        --outfile spades_norm${norm}_k${kmer_size}_\${subset_id}.fasta \
        --mode rna \
        --kmer_size ${kmer_size} \
        --strandedness ${params.strandedness}
    
    mv fofn.txt fofn_norm${norm}_k${kmer_size}_\${subset_id}.txt
    
    cp .command.log logs/slurm_norm${norm}_k${kmer_size}_\${subset_id}.log
    """
}
