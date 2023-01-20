process EVIGENE {
    tag "Merge all assemblies into a single FASTA file"
    publishDir "${params.outdir}/evigene", mode: 'copy'

    input:
    path concat_assembly

    output:
    path "out/evigene.fasta", emit: assembly_alltrans
    path "out/evigene_primarytrans.fasta", emit: assembly_1trans
    path "out/okayset", emit: okayset
    path "out/logs/slurm*log", emit: log
    
    script:
    """
    evigene.sh \
        --infile ${concat_assembly} \
        --outfile out/evigene.fasta

    cp .command.log out/logs/slurm.log
    """
}
