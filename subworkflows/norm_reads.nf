// Subworkflow to create transcriptome assemblies with normalized reads

// Include modules
include { ORNA_ALL; ORNA_PAIR } from '../modules/orna'
include { FOFN } from '../modules/fofn'

// Define the subworkflow
workflow NORM_READS {
    take:
        Reads
    
    main:
        // Read normalization with ORNA
        ch_orna_pair = ORNA_PAIR(Reads).fq
        ch_reads_fofn = FOFN(ch_orna_pair.collect().flatten().collect(), "fastq.gz")
        ch_reads = ORNA_ALL(ch_reads_fofn).fq
    
    emit:
        ch_reads
}
