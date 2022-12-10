process EVIGENE {
    tag "Merge all assemblies into a single FASTA file"
    publishDir "${params.outdir}/evigene", mode: 'copy'

    input:
    path concat_assembly

    output:
    path "out/final/evigene_all.fasta", emit: assembly_alltrans
    path "out/final/evigene_primarytrans.fasta", emit: assembly_1trans
    path "out/okayset", emit: okayset
    path "out/logs/slurm*log", emit: log
    
    script:
    """
    evigene.sh \
        --infile ${concat_assembly} \
        --outdir out

    cp .command.log out/logs/slurm.log
    """
}
