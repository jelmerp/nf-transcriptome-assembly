process KALLISTO_INDEX {
    tag "Index the final transcriptome assembly with Kallisto"
    publishDir "${params.outdir}/kallisto_index", mode: 'copy'

    input:
    path merged_assembly

    output:
    path "assembly.idx", emit: index
    path "logs/slurm*log", emit: log
    
    script:
    """
    kallisto_index.sh \
        --infile ${merged_assembly} \
        --outfile assembly.idx

    cp .command.log logs/slurm.log
    """
}
