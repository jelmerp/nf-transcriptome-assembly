// Normalize all samples together
process ORNA_ALL {
    publishDir "${params.outdir}/orna", mode: 'copy'

    input:
    path fofn

    output:
    path "orna_out/*fastq.gz", emit: fq
    path "orna_out/logs/slurm*log", emit: log
    path "orna_out/logs/version.txt", emit: version
    
    script:
    """
    orna.sh \
        --fofn ${fofn} \
        --outdir orna_out
    
    cp .command.log orna_out/logs/slurm.log
    """
}

// Normalize one sample
process ORNA_PAIR {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fq_pair)

    output:
    path "orna_out/*fastq.gz", emit: fq
    path "orna_out/logs/slurm*log", emit: log
    path "orna_out/logs/version.txt", emit: version
    //Need to go into outdir 'orna_out' or it will try to overwrite the input files

    script:
    """
    orna.sh \
        --R1 ${fq_pair[0]} \
        --outdir orna_out
    
    cp .command.log orna_out/logs/slurm-${sample_id}.log
    """
}
