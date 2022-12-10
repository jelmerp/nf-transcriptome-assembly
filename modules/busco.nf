
process BUSCO {
    tag "Evaluate an assembly with Busco"
    publishDir "${params.outdir}/busco", mode: 'copy'

    input:
    path assembly
    val busco_db

    output:
    path "*", emit: out
    path "logs/slurm*log", emit: log
    
    script:
    """
    busco.sh \
        --infile ${assembly} \
        --db ${busco_db} \
        --outdir .

    cp .command.log logs/slurm.log
    """
}
