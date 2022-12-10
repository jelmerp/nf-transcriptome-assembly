// Subworkflow to process raw RNAseq reads prior to transcriptome assembly

// Parameters assumed to be passed from main.nf:
// - params.nfiles_rcorr
// - params.kraken_db_url   (Can be 'false')
// - params.kraken_db_dir   (Can be 'false')

// Include modules
include { MULTIQC as MULTIQC_TRIM; MULTIQC as MULTIQC_PRE } from '../modules/multiqc'
include { TRIMGALORE } from '../modules/trimgalore'
include { FOFN_FOR_RCORRECTOR } from '../modules/fofn_for_rcorrector'
include { RCORRECTOR } from '../modules/rcorrector'
include { RCORRFILTER } from '../modules/rcorrfilter'
include { SORTMERNA } from '../modules/sortmerna'
include { GET_KRAKEN_DB } from '../modules/get_kraken_db'
include { KRAKEN; JOIN_KRAKEN } from '../modules/kraken'
include { KRONA } from '../modules/krona'
include { ORNA; JOIN_ORNA } from '../modules/orna'

// Define the subworkflow
workflow PROCESS_READS {
    take:
        Reads

    main:
        // Trimming with TrimGalore
        ch_trim = TRIMGALORE(Reads)
        MULTIQC_TRIM(ch_trim.fqc_zip.mix(ch_trim.trim_report).collect(), "multiqc_trimmed.html")

        // Read correction with Rocorrector
        ch_trim_all = CAT_FILENAMES(ch_trim.fq_trimmed.collect().flatten().collect())
        ch_rcorr = RCORRECTOR(ch_trim_all.splitText(by: params.nfiles_rcorr)) // Run with nfiles_rcorr files per time

        ch_rcorr = ch_rcorr.fq             // Get reads back into by-sample tuple format
            .flatten()
            .map { file -> tuple(file.simpleName - ~/_R[12].*/, file) }
            .groupTuple(by: 0, size: 2)
        ch_rcorrfilter = RCORRFILTER(ch_rcorr)
        
        // Remove rRNA with SortMeRNA
        ch_sortmerna = SORTMERNA(ch_rcorrfilter.fq)
        
        // Remove contaminants with Kraken
        if (params.kraken_db_url != false) {
            ch_kraken_db = GET_KRAKEN_DB(params.kraken_db_url)
        } else {
            ch_kraken_db = Channel.fromPath(params.kraken_db_dir)
        }
        ch_kraken = KRAKEN(ch_sortmerna.fq_unmapped, ch_kraken_db)
        KRONA(ch_kraken.main)
        MULTIQC_PRE(ch_kraken.report.mix(ch_sortmerna.log).collect(), "multiqc_preprocess")

        // Read normalization with ORNA
        ch_orna = ORNA(ch_kraken.fq)

        // Get channels with all preprocessed FASTQ files in one dir
        ch_normreads = JOIN_ORNA(ch_orna.fq.collect())
        ch_nonormreads = JOIN_KRAKEN(ch_kraken.fq_list.collect())
    
    emit:
        reads_norm = ch_normreads
        reads_nonorm = ch_nonormreads
}
