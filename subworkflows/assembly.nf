// Subworkflow to assemble a transcriptome

// Include modules
include { INDEX_GENOME } from '../modules/index_genome'
include { MAP2GENOME } from '../modules/map2genome'
include { MERGE_BAM } from '../modules/merge_bam'
include { TRINITY as TRINITY_NORM; TRINITY as TRINITY_NONORM } from '../modules/trinity'
include { TRINITY_GUIDED } from '../modules/trinity_guided'
include { TRANSABYSS as TRANSABYSS_NORM; TRANSABYSS as TRANSABYSS_NONORM } from '../modules/transabyss'
include { SPADES as SPADES_NORM; SPADES as SPADES_NONORM } from '../modules/spades'
include { CONCAT_ASSEMBLIES } from '../modules/concat_assemblies'
include { EVIGENE } from '../modules/evigene'
include { FOFN as FOFN_BAM; FOFN as FOFN_NONORM } from '../modules/fofn'

// Process params
nfiles = params.nfiles_assembly
genome_fasta = params.genome_fasta
genome_index = params.genome_index
genome_bam = params.genome_bam
assemble_norm = params.skip_assemble_norm == true ? false : true
assemble_nonorm = params.skip_assemble_nonorm == true ? false : true
assemble_trinity = params.skip_assemble_trinity == true ? false : true
assemble_spades = params.skip_assemble_spades == true ? false : true
assemble_abyss = params.skip_assemble_abyss == true ? false : true

// Define the subworkflow
workflow ASSEMBLY {
    take:
        Reads_norm
        Reads_nonorm
        K_abyss
        K_spades
    
    main:
        // Split FASTQ files into groups
        ch_reads_nonorm = FOFN_NONORM(Reads_nonorm.collect().flatten().collect(), "fastq.gz")
        ch_reads_nonorm = ch_reads_nonorm.splitText(by: nfiles)

        // Reference-guided Trinity
        if (assemble_trinity == true && (genome_fasta != false || genome_index != false)) {
            // Get or create index
            if (genome_index == false) {
                ch_genome = channel.fromPath(genome_fasta)
                ch_index = INDEX_GENOME(ch_genome).index
            } else {
                ch_index = channel.value(genome_index) // NOTE: Channel.fromPath() doesn't work here - won't recycle
            }
            // Get or create BAM
            if (genome_bam == false) {
                ch_bam = MAP2GENOME(Reads_nonorm, ch_index)
                ch_bam_fofn = FOFN_BAM(ch_bam.bam.collect(), "bam")
                ch_bam_grouped = ch_bam_fofn.splitText(by: nfiles)
                ch_bam_merged = MERGE_BAM(ch_bam_grouped).bam
            } else {
                ch_bam_merged = channel.fromPath(genome_bam)
            }
            // Reference-guided assembly
            ch_trinity_gg = TRINITY_GUIDED(ch_bam_merged).assembly
        } else {
            ch_trinity_gg = channel.empty()
        }

        // De novo Trinity
        ch_trinity_norm = channel.empty()
        ch_trinity_nonorm = channel.empty()
        if (assemble_trinity == true) {
            if (assemble_norm == true) ch_trinity_norm = TRINITY_NORM(Reads_norm, "true").assembly
            if (assemble_nonorm == true) ch_trinity_nonorm = TRINITY_NONORM(ch_reads_nonorm, "false").assembly
        }

        // Trans-abyss
        ch_abyss_norm = channel.empty()
        ch_abyss_nonorm = channel.empty()
        if (assemble_abyss == true) {
            if (assemble_norm == true) ch_abyss_norm = TRANSABYSS_NORM(Reads_norm, K_abyss, "true").assembly
            if (assemble_nonorm == true) ch_abyss_nonorm = TRANSABYSS_NONORM(ch_reads_nonorm, K_abyss, "false").assembly
        }

        // SPAdes
        ch_spades_norm = channel.empty()
        ch_spades_nonorm = channel.empty()
        if (assemble_spades == true) {
            if (assemble_norm == true) ch_spades_norm = SPADES_NORM(Reads_norm, K_spades, "true").assembly
            if (assemble_nonorm == true) ch_spades_nonorm = SPADES_NONORM(ch_reads_nonorm, K_spades, "false").assembly
        }

        // Merge assemblies
        ch_assemblies_all = ch_trinity_gg
            .mix(ch_trinity_norm, ch_trinity_nonorm, ch_abyss_norm, ch_abyss_nonorm, ch_spades_norm, ch_spades_nonorm)
            .collect()

        ch_assemblies_concat = CONCAT_ASSEMBLIES(ch_assemblies_all)
        ch_evigene = EVIGENE(ch_assemblies_concat.assembly)

    emit:
        assembly_1trans = ch_evigene.assembly_1trans
        assembly_alltrans = ch_evigene.assembly_alltrans
}
