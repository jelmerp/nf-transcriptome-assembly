// Subworkflow to assemble a transcriptome

// Parameters assumed to be passed from main.nf:
// - params.genome_fasta
// - params.genome_index
// - k_abyss
// - k_spades

// Include modules
include { INDEX_GENOME } from '../modules/index_genome'
include { MAP2GENOME } from '../modules/map2genome'
include { MERGE_BAM } from '../modules/merge_bam'
include { TRINITY } from '../modules/trinity'
include { TRINITY_GUIDED } from '../modules/trinity_guided'
include { TRANSABYSS as TRANSABYSS_NORM; TRANSABYSS as TRANSABYSS_NONORM } from '../modules/transabyss'
include { SPADES as SPADES_NORM; SPADES as SPADES_NONORM } from '../modules/spades'
include { CONCAT_ASSEMBLIES } from '../modules/concat_assemblies'
include { EVIGENE } from '../modules/evigene'

// Process parameters
if (params.genome_bam == true) {
    ch_bam_merged = Channel.fromPath(params.genome_bam)
    skip_map2genome = true
}

// Hardcoded parameters
norms = [ 'true', 'false']    // Run assemblies with normalized and non-normalized data

// Define the subworkflow
workflow ASSEMBLY {
    take:
        Reads_norm
        Reads_nonorm
        K_abyss
        K_spades
    
    main:
        // Reference-guided Trinity
        if (params.genome_fasta != false || params.genome_index != false ) {
            if (params.genome_index == false) {
                ch_genome = Channel.fromPath(params.genome_fasta)
                ch_index = INDEX_GENOME(ch_genome).index
            } else {
                ch_index = channel.value(params.genome_index) // NOTE: Channel.fromPath() doesn't work here - won't recycle
            }
            if (skip_map2genome == false) {
                ch_bam = MAP2GENOME(Reads_nonorm, ch_index)
                ch_bam_merged = MERGE_BAM(ch_bam.bam.collect())
            }
            ch_trinity_gg = TRINITY_GUIDED(ch_bam_merged.bam)
        } else {
            ch_trinity_gg = channel.empty()
        }

        // De novo Trinity
        ch_trinity = TRINITY(Reads_nonorm, norms).assembly

        // Trans-abyss
        ch_abyss_norm = TRANSABYSS_NORM(Reads_norm, K_abyss, "true").assembly
        ch_abyss_nonorm = TRANSABYSS_NONORM(Reads_nonorm, K_abyss, "false").assembly

        // SPAdes
        ch_spades_norm = SPADES_NORM(Reads_norm, K_spades, "true").assembly
        ch_spades_nonorm = SPADES_NONORM(Reads_nonorm, K_spades, "false").assembly

        // Merge assemblies
        ch_assemblies_all = ch_trinity_gg
            .mix(ch_trinity, ch_abyss_norm, ch_abyss_nonorm, ch_spades_norm, ch_spades_nonorm)
            .collect()
        ch_assemblies_concat = CONCAT_ASSEMBLIES(ch_assemblies_all)
        ch_evigene = EVIGENE(ch_assemblies_concat.assembly)

    emit:
        assembly_1trans = ch_evigene.assembly_1trans
        assembly_alltrans = ch_evigene.assembly_alltrans
}
