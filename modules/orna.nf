process ORNA {
    tag "ORNA read normalization for ${sample_id}"
    publishDir "${params.outdir}/orna", mode: 'copy'

    input:
    path fq_pair

    output:
    path "orna_out/*fastq.gz", emit: fq
    path "orna_out/logs/slurm*log", emit: log
    path "orna_out/logs/version.txt", emit: version
    
    script:
    """
    orna.sh \
        --R1 ${fq_pair[0]} \
        --outdir orna_out
    
    cp .command.log orna_out/logs/slurm-${sample_id}.log
    """
}

//TODO can probably delete this
process JOIN_ORNA {
    tag "Combining ORNA output for assembly processes"
    label "local_process"

    input:
    path orna_fastqs

    output:
    path "orna_fastqs_all"

    script:
    """
    mkdir orna_fastqs_all
    mv -v $orna_fastqs orna_fastqs_all
    """
}
