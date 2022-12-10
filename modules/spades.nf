process SPADES {
    tag "Assembly with SPADES - kmer $kmer_size - normalization $norm"
    publishDir "${params.outdir}/spades", mode: 'copy'

    input:
    path dir_with_all_fqs
    each kmer_size
    val norm

    output:
    path "spades_norm*.fasta", emit: assembly
    path "logs/slurm*log", emit: log
    path "logs/version.txt", emit: version
    
    script:
    """
    spades.sh \
        --indir ${dir_with_all_fqs} \
        --outdir . \
        --mode rna \
        --kmer_size ${kmer_size} \
        --strandedness ${params.strandedness}
    
    mv transcripts.fasta spades_norm${norm}_k${kmer_size}.fasta

    cp .command.log logs/slurm.log
    """
}
