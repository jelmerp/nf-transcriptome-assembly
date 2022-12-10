process TXIMPORT {
    tag "Import Kallisto counts into R with tximport"
    publishDir "${params.outdir}/kallisto", mode: 'copy'

    input:
    path kallisto_files
    path tx2gene

    output:
    path "tximport/tximport_object.rds", emit: tximport
    path "tximport/logs/slurm*log", emit: log
    
    script:
    """
    kallisto_all_dir=kallisto_all

    mkdir -p \${kallisto_all_dir}
    mv -v $kallisto_files \${kallisto_all_dir}

    tximport.R \
        --indir \${kallisto_all_dir} \
        --tx2gene ${tx2gene} \
        --outdir tximport

    cp .command.log tximport/logs/slurm.log
    """
}
