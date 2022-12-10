process RNAQUAST {
    tag "Evaluate an assembly with RNAQuast"
    publishDir "${params.outdir}/rnaquast", mode: 'copy'

    input:
    path assembly

    output:
    path "*report*", emit: out
    path "logs/slurm*log", emit: log
    
    script:
    """
    rnaquast.sh \
        --assemblies ${assembly} \
        --outdir . \
        --strandedness ${params.strandedness}

    cp .command.log logs/slurm.log
    """
}
