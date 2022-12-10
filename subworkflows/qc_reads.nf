// Subworkflow to QC reads in FASTQ files

// Include modules
include { FASTQC } from '../modules/fastqc'
include { MULTIQC } from '../modules/multiqc'

// Define the subworkflow
workflow QC_READS {
    take:
        Reads

    main:
        ch_fastqc = FASTQC(Reads)
        MULTIQC(ch_fastqc.zip.collect(), "multiqc_raw")
}
