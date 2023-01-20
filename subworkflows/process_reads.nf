// Subworkflow to process raw RNAseq reads prior to transcriptome assembly

// Include modules
include { MULTIQC as MULTIQC_TRIM; MULTIQC as MULTIQC_PRE } from '../modules/multiqc'
include { TRIMGALORE } from '../modules/trimgalore'
include { FOFN } from '../modules/fofn'
include { RCORRECTOR } from '../modules/rcorrector'
include { RCORRFILTER } from '../modules/rcorrfilter'
include { SORTMERNA } from '../modules/sortmerna'
include { GET_KRAKEN_DB } from '../modules/get_kraken_db'
include { KRAKEN } from '../modules/kraken'
include { KRONA } from '../modules/krona'
include { MV_IN_1DIR } from '../modules/mv_in_1dir'

// Process params
nfiles = params.nfiles_rcorr
kraken_db_url = params.kraken_db_url
kraken_db_dir = params.kraken_db_dir

// Define the subworkflow
workflow PROCESS_READS {
    take:
        Reads

    main:
        // Trimming with TrimGalore
        ch_trim = TRIMGALORE(Reads)
        MULTIQC_TRIM(ch_trim.fqc_zip.mix(ch_trim.trim_report).collect(), "multiqc_trimmed.html")

        // Read correction with Rcorrector
        ch_trim_all = FOFN(ch_trim.fq_trimmed.collect().flatten().collect(), "fastq.gz")
        ch_rcorr = RCORRECTOR(ch_trim_all.splitText(by: nfiles)) // Run with nfiles_rcorr files per time

        ch_rcorr = ch_rcorr.fq             // Get reads back into by-sample tuple format
            .flatten()
            .map { file -> tuple(file.simpleName - ~/_R[12].*/, file) }
            .groupTuple(by: 0, size: 2)
        ch_rcorrfilter = RCORRFILTER(ch_rcorr)
        
        // Remove rRNA with SortMeRNA
        ch_sortmerna = SORTMERNA(ch_rcorrfilter.fq)
        
        // Remove contaminants with Kraken
        if (kraken_db_url != false) {
            ch_kraken_db = GET_KRAKEN_DB(kraken_db_url)
        } else {
            ch_kraken_db = Channel.fromPath(kraken_db_dir)
        }
        ch_kraken = KRAKEN(ch_sortmerna.fq_unmapped, ch_kraken_db)
        
        // Plot Kraken output with Krona
        KRONA(ch_kraken.main)
        
        // Run MultiQC on the output of all previous tools
        MULTIQC_PRE(ch_kraken.report.mix(ch_sortmerna.log).collect(), "multiqc_preprocess")
    
    emit:
        reads = ch_kraken.fq
}
