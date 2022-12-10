process CONCAT_ASSEMBLIES {
    tag "Concatenate all assemblies into a single FASTA file"
    publishDir "${params.outdir}/concat_assembly", mode: 'copy'

    input:
    path assemblies

    output:
    path "concat_assembly.fasta", emit: assembly
    path "logs/slurm*log", emit: log
    
    script:
    """
    concat_assemblies.sh \
        --outfile concat_assembly.fasta \
        ${assemblies}

    cp .command.log logs/slurm.log
    """
}
