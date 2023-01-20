// Subworkflow to create transcriptome assemblies with non-normalized reads

// Include modules
include { TRINITY } from '../modules/trinity'
include { TRANSABYSS } from '../modules/transabyss'
include { SPADES } from '../modules/spades'
include { FOFN as FOFN_READS } from '../modules/fofn'

// Process params
nfiles = params.nfiles_assembly
run_trinity = params.skip_trinity == true ? false : true
run_spades = params.skip_spades == true ? false : true
run_abyss = params.skip_abyss == true ? false : true

// Default empty channels
ch_trinity = channel.empty()
ch_abyss = channel.empty()
ch_spades = channel.empty()

// Define the subworkflow
workflow ASSEMBLE_NONORM {
    take:
        Reads
        K_abyss
        K_spades
    
    main:
        // Split non-FASTQ files into groups
        ch_read_fofn = FOFN_READS(Reads.collect(), "fastq.gz") //.flatten().collect()
        ch_reads = ch_read_fofn.splitText(by: nfiles)

        // Trans-abyss
        //if (run_abyss) ch_abyss = TRANSABYSS(ch_reads, K_abyss, "false").assembly

        // SPAdes
        //if (run_spades) ch_spades = SPADES(ch_reads, K_spades, "false").assembly

        // Trinity
        if (run_trinity) ch_trinity = TRINITY(ch_reads, "false").assembly

        // Combine assemblies into a single channel
        ch_assemblies_all = ch_trinity.mix(ch_abyss, ch_spades).collect()

    emit:
        ch_assemblies_all
}
