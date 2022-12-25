// Subworkflow to process raw RNAseq reads prior to transcriptome assembly

// Include modules
include { MULTIQC as MULTIQC_TRIM; MULTIQC as MULTIQC_PRE } from '../modules/multiqc'
include { TRIMGALORE } from '../modules/trimgalore'
include { FOFN_CREATE } from '../modules/fofn_create'
include { RCORRECTOR } from '../modules/rcorrector'
include { RCORRFILTER } from '../modules/rcorrfilter'
include { SORTMERNA } from '../modules/sortmerna'
include { GET_KRAKEN_DB } from '../modules/get_kraken_db'
include { KRAKEN; JOIN_KRAKEN } from '../modules/kraken'
include { KRONA } from '../modules/krona'
include { ORNA; JOIN_ORNA } from '../modules/orna'
include { CONCAT_FASTQS } from '../modules/concat_fastqs'

// Process params
assemble_norm = params.skip_assemble_norm == true ? false : true
assemble_nonorm = params.skip_assemble_nonorm == true ? false : true
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
        ch_trim_all = FOFN_CREATE(ch_trim.fq_trimmed.collect().flatten().collect())
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
        
        KRONA(ch_kraken.main)
        
        MULTIQC_PRE(ch_kraken.report.mix(ch_sortmerna.log).collect(), "multiqc_preprocess")

        // Read normalization with ORNA
        if (assemble_norm == true) {
            ch_concat = CONCAT_FASTQS(ch_kraken.fq.collect())
            ch_orna = ORNA(ch_concat.fq)
        }

        // Get channels with all preprocessed FASTQ files in one dir
        ch_normreads = channel.empty()
        ch_nonormreads = channel.empty()
        if (assemble_norm == true) ch_normreads = JOIN_ORNA(ch_orna.fq.collect())
        if (assemble_nonorm == true) ch_nonormreads = JOIN_KRAKEN(ch_kraken.fq_list.collect())
    
    emit:
        reads_norm = ch_normreads
        reads_nonorm = ch_nonormreads
}
