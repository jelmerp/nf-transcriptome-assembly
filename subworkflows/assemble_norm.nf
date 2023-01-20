// Subworkflow to create transcriptome assemblies with normalized reads

// Include modules
include { SPADES as SPADES_NORM } from '../modules/spades'
include { TRANSABYSS as TRANSABYSS_NORM } from '../modules/transabyss'
include { TRINITY as TRINITY_NORM } from '../modules/trinity'
include { INDEX_GENOME } from '../modules/index_genome'
include { MAP2GENOME } from '../modules/map2genome'
include { MERGE_BAM } from '../modules/merge_bam'
include { TRINITY_GUIDED } from '../modules/trinity_guided'
include { FOFN } from '../modules/fofn'

// Process params
genome_fasta = params.genome_fasta
genome_index = params.genome_index
genome_bam = params.genome_bam
nfiles = params.nfiles_assembly

// Determine what to run
run_trinity = params.skip_trinity == true ? false : true
run_spades = params.skip_spades == true ? false : true
run_abyss = params.skip_abyss == true ? false : true
if (run_trinity == true && (genome_fasta != false || genome_index != false)) {
    run_trinity_gg = true
} else {
    run_trinity_gg = false
}

// Default empty channels
ch_abyss = channel.empty()
ch_spades = channel.empty()
ch_trinity = channel.empty()
ch_trinity_gg = channel.empty()

// Define the subworkflow
workflow ASSEMBLE_NORM {
    take:
        Reads
        K_abyss
        K_spades
    
    main:
        // Get reads in right format
        ch_read_fofn = FOFN(Reads.collect(), "fastq.gz")
        ch_reads = ch_read_fofn.splitText(by: nfiles)

        // Trans-abyss
        if (run_abyss) ch_abyss = TRANSABYSS_NORM(ch_reads, K_abyss, "true").assembly

        // SPAdes
        if (run_spades) ch_spades = SPADES_NORM(ch_reads, K_spades, "true").assembly

        // De novo Trinity
        if (run_trinity) ch_trinity = TRINITY_NORM(ch_reads, "true").assembly

        // Reference-guided Trinity
        if (run_trinity_gg == true) {
            
            // Get or create genome index
            if (genome_index == false) {
                ch_genome = channel.fromPath(genome_fasta, checkIfExists: true)
                ch_index = INDEX_GENOME(ch_genome)
                ch_index = ch_index.index.collect() // Needed to make Nextflow recycle the index for each sample
            } else {
                ch_index = channel.fromPath(genome_index, checkIfExists: true)
                ch_index = ch_index.collect()  // Needed to make Nextflow recycle the index for each sample
            }

            // Get or create BAM
            if (genome_bam == false) {
                ch_bam = MAP2GENOME(ch_reads, ch_index)
                ch_bam_merged = MERGE_BAM(ch_bam.bam.collect()).bam
                //ch_bam_fofn = FOFN(ch_bam.bam.collect(), "bam")
                //ch_bam_grouped = ch_bam_fofn.splitText(by: nfiles)
                //ch_bam_merged = MERGE_BAM(ch_bam_grouped).bam
            } else {
                ch_bam_merged = channel.fromPath(genome_bam, checkIfExists: true)
            }
            
            // Reference-guided assembly
            ch_trinity_gg = TRINITY_GUIDED(ch_bam_merged).assembly
        }

        // Combine assemblies into a single channel
        ch_assemblies_all = ch_trinity_gg.mix(ch_trinity, ch_abyss, ch_spades).collect()

    emit:
        ch_assemblies_all
}
