
process KALLISTO_QUANT {
    tag "Quantify read counts with Kallisto for ${sample_id}"
    publishDir "${params.outdir}/kallisto", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)
    path assembly_index

    output:
    path "*abundance.h5", emit: counts_h5
    path "*abundance.tsv", emit: counts_tsv
    path "${sample_id}/logs/slurm*log", emit: log
    
    script:
    """
    kallisto_quant.sh \
        --R1 ${fq_pair[0]} \
        --ref-index ${assembly_index} \
        --strandedness ${params.strandedness} \
        --outdir .

    cp .command.log ${sample_id}/logs/slurm.log
    """
}
