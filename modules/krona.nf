process KRONA {
    tag "Krona plots of Kraken results for ${sample_id}"
    publishDir "${params.outdir}/krona", mode: 'copy'

    input:
    tuple val(sample_id), path(kraken_res)

    output:
    path "*html"
    path "logs/slurm*log"
    
    script:
    """
    krona.sh \
        --infile ${kraken_res} \
        --outfile ./krona_${sample_id}.html
    
    cp .command.log logs/slurm-${sample_id}.log
    """
}
