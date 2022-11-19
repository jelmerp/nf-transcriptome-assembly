#!/usr/bin/env nextflow

// Nextflow transcriptome assembly pipeline

def helpMessage() {
    log.info"""
    ============================================================================
            T R A N S C R I P T O M E   A S S E M B L Y   P I P E L I N E
    ============================================================================
    REQUIRED OPTIONS:
        --reads             <str>   Single-quoted stirng with path to input dir + glob to match FASTQ files
                                    E.g. 'data/fastq/*_R{1,2}.fastq.gz'
        --busco_db          <str>   Busco Database name
    
    OTHER OPTIONS:
        --outdir            <dir>   Final output dir for workflow results               [default: 'results/nf_tram']
        --subset_fastq      <int>   Subset the FASTQ files to <int> reads               [default: no subsetting]
        --ref_fasta         <file>  Reference FASTA file for Trinity genome-guided assembly [default: unset]
        --trim_nextseq              Use TrimGalore/Cutadapt's 'nextseq' option for poly-G trimming [default: don't use]
        --nfiles_rcorr      <int>   Number of files to run rcorrector with at a time    [default: 20 (=10 PE samples)]
        --strandedness      <str>   'reverse', 'forward', or 'unstranded'               [default: 'reverse']
        --min_contig_length <int>   Minimum contig length (TransAbyss and Trinity)      [default: 300]
        --k_vals_transabyss <str>   Comma-separated list of kmer-values for TransAbyss
        --k_vals_spades     <str>   Comma-separated list of kmer-values for SPAdes
        --kraken_db         <dir>   Path to a Kraken database                           [default: '/fs/project/PAS0471/jelmer/refdata/kraken/std-plus-fungi']
        --entap_config      <file>  EnTAP Config file
    UTILITY OPTIONS
        --help                      Print this help message and exit
    """.stripIndent()
}

//  print help if required
if (params.help) {
    helpMessage()
    exit 0
}

// Process parameters
if (!params.reads) { exit 1, '\n============\nERROR: Input reads not specified! Use "--reads" to do so\n============\n' }
if (!params.outdir) { exit 1, '\n============\nERROR: Output dir not specified! Use "--outdir" to do so\n============\n' }

def checkPathParamList = [ params.reads ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

n_reads = params.subset_fastq
k_vals_abyss = params.k_vals_abyss
k_vals_spades = params.k_vals_spades

// Hardcoded parameters
norms = [ 'true', 'false']    // Run assemblies with normalized and non-normalized data

// Run info
log.info ""
log.info "======================================================================"
log.info "    T R A N S C R I P T O M E   A S S E M B L Y   P I P E L I N E"
log.info "======================================================================"
log.info "Parameters:"
params.each { k, v -> if (v) { println "${k}: ${v}" } }
log.info "======================================================================"
log.info ""

// Include module
include { SUBSET_FASTQ } from './modules/all_mods'
include { FASTQC } from './modules/all_mods'
include { MULTIQC_RAW; MULTIQC_TRIM; MULTIQC_PRE } from './modules/all_mods'
include { TRIMGALORE; CAT_FILENAMES } from './modules/all_mods'
include { RCORRECTOR; RCORRFILTER } from './modules/all_mods'
include { SORTMERNA } from './modules/all_mods'
include { KRAKEN; KRONA } from './modules/all_mods'
include { ORNA; JOIN_ORNA; JOIN_KRAKEN } from './modules/all_mods'
include { INDEX_GENOME; MAP2GENOME } from './modules/all_mods'
include { TRINITY ; TRINITY_GUIDED } from './modules/all_mods'
include { TRANSABYSS as TRANSABYSS_NORM; TRANSABYSS as TRANSABYSS_NONORM } from './modules/all_mods'
include { SPADES as SPADES_NORM; SPADES as SPADES_NONORM } from './modules/all_mods'
//include { CONCAT_ASSEMBLIES; EVIGENE } from './modules/all_mods'
//include { TRINITY_STATS; BUSCO; RNAQUAST; DETONATE } from './modules/all_mods'
//include { ENTAP_CONFIG; ENTAP } from './modules/all_mods'
//include { KALLISTO_INDEX; KALLISTO } from './modules/all_mods'

// Workflow
workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    if (params.ref_fasta != false) {
        ref_ch = Channel.fromPath(params.ref_fasta)
    }

    println "Number of input samples:"
    reads_ch.count().view()
    
    // ======================================================================= //
    //                            PRE-PROCESSING
    // ======================================================================= //
    if (params.subset_fastq != false) {
        finalreads_ch = SUBSET_FASTQ(reads_ch, n_reads).fq
    } else {
        finalreads_ch = reads_ch
    }

    // Initial QC
    fqc_ch = FASTQC(finalreads_ch)
    MULTIQC_RAW(fqc_ch.zip.collect())

    // Trimming
    trim_ch = TRIMGALORE(finalreads_ch)
    MULTIQC_TRIM(trim_ch.fqc_zip.mix(trim_ch.trim_report).collect())

    // Read correction
    trim_all_ch = CAT_FILENAMES(trim_ch.fq_trimmed.collect().flatten().collect())
    rcorr_ch = RCORRECTOR(trim_all_ch.splitText(by: params.nfiles_rcorr)) // Run with nfiles_rcorr files per time

    rcorr_ch = rcorr_ch.fq             // Get reads back into by-sample tuple format
        .flatten()
        .map { file -> tuple(file.simpleName - ~/_R[12].*/, file) }
        .groupTuple(by: 0, size: 2)
    rcorrfilter_ch = RCORRFILTER(rcorr_ch)
    
    // Remove rRNA
    sortmerna_ch = SORTMERNA(rcorrfilter_ch.fq)
    
    // Remove contaminants
    kraken_ch = KRAKEN(sortmerna_ch.fq_unmapped)
    KRONA(kraken_ch.main)
    MULTIQC_PRE(kraken_ch.report.mix(sortmerna_ch.log).collect())

    // Read normalization
    orna_ch = ORNA(kraken_ch.fq)
    
    // ======================================================================= //
    //                              ASSEMBLY
    // ======================================================================= //
    // Get channels with all preprocessed FASTQ files in one dir
    normreads_ch = JOIN_ORNA(orna_ch.fq.collect())
    nonormreads_ch = JOIN_KRAKEN(kraken_ch.fq_list.collect())

    // Reference-guided Trinity
    if (params.ref_fasta != false) {
        index_ch = INDEX_GENOME(ref_ch)
        bam_ch = MAP2GENOME(ref_ch, kraken_ch.fq)
        trinity_gg_ch = TRINITY_GUIDED(bam_ch.bam)
    }

    // De novo Trinity
    trinity_nonorm_ch = TRINITY(nonormreads_ch, norms)

    // Trans-abyss
    transabyss_norm_ch = TRANSABYSS_NORM(normreads_ch, k_vals_abyss, "true")
    transabyss_nonorm_ch = TRANSABYSS_NONORM(nonormreads_ch, k_vals_abyss, "false")

    // SPAdes
    spades_norm_ch = SPADES_NORM(normreads_ch, k_vals_spades, "true")
    spades_nonorm_ch = SPADES_NONORM(nonormreads_ch, k_vals_spades, "true")

    // Merge assemblies
    // concat_asm_ch = CONCAT_ASSEMBLIES()
    // evigene_ch = EVIGENE(concat_asm_ch)

    // QC assemblies
    //TRINITY_STATS()
    //BUSCO()
    //RNAQUAST()
    //DETONATE()

    // Annotation
    // ENTAP_CONFIG(params.entap_config)
    // ENTAP(evigene_ch.assembly_primarytrans, params.entap_config)

    // Quantification
    // kallisto_idx_ch = KALLISTO_INDEX(evigene_ch)
    // KALLISTO(kraken_ch, kallisto_idx_ch)
}

// Report completion/failure of workflow
workflow.onComplete {
    println ( workflow.success ? """
        ========================================================================
                            PIPELINE SUCCESSFULLY COMPLETED
        ========================================================================
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        ========================================================================
                                    PIPELINE FAILED
        ========================================================================
        """
        .stripIndent()
    )
}
