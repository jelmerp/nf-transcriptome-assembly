process TRANSABYSS {
    tag "Assembly with Trans-ABySS - kmer $kmer_size - normalization $norm"
    publishDir "${params.outdir}/transabyss", mode: 'copy'

    input:
    path dir_with_all_fqs
    each kmer_size
    val norm

    output:
    path "transabyss_norm*.fasta", emit: assembly
    path "logs/slurm*log", emit: log
    path "logs/version.txt", emit: version
    
    script:
    """
    transabyss.sh \
        --indir ${dir_with_all_fqs} \
        --outdir . \
        --id transabyss_norm${norm}_k${kmer_size} \
        --min_contig_length ${params.min_contig_length} \
        --kmer_size ${kmer_size} \
        --strandedness ${params.strandedness}
    
    mv -v transabyss_norm${norm}_k${kmer_size}-final.fa transabyss_norm${norm}_k${kmer_size}.fasta

    cp .command.log logs/slurm.log
    """
}
