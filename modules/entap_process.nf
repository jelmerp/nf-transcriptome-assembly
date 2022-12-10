process ENTAP_PROCESS {
    tag "Post-processing of EnTap output"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path entap_out
    path assembly_1trans
    path assembly_alltrans

    output:
    path "entap_ed", emit: out
    path "entap_ed/assembly.fasta", emit: assembly
    path "entap_ed/gene_trans_map.tsv", emit: tx2gene
    path "logs/slurm*log", emit: log
    
    script:
    """
    entap_process.sh \
        --entap-dir ${entap_out} \
        --in-1trans ${assembly_1trans} \
        --in-alltrans ${assembly_alltrans} \
        --outdir entap_ed

    cp .command.log entap_ed/logs/slurm.log
    """
}
