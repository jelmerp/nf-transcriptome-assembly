#!/usr/bin/env nextflow

// Nextflow transcriptome assembly pipeline

//TODO
// - TransABYSS: copy input to $TMPDIR for faster reading
// - Split mapping of normalized reads - split FASTQ files
// - All modules: process 'tags' only for processes that run multiple times + should be succinct
// - Don't use local Conda envs

// Help function
def helpMessage() {
    log.info"""
    ===========================================================================================================================
                            T R A N S C R I P T O M E   A S S E M B L Y   P I P E L I N E
    ===========================================================================================================================
    REQUIRED OPTIONS:
        --reads                 <str>   Single-quoted string with path to input dir + glob to match FASTQ files
                                            E.g. 'data/fastq/*_R{1,2}.fastq.gz'
        --busco_db              <str>   Busco database name (see https://busco.ezlab.org/list_of_lineages.html)
        --entap_taxon           <file>  Taxon name for EnTAP, using format 'homo_sapiens' (lower case and underscores instead of spaces)
        --entap_contam          <file>  Comma-separated list of contaminant taxa for EnTAP (e.g., 'viruses,bacteria')
                                            Taxon names should be according to NCBI taxonomy (https://www.ncbi.nlm.nih.gov/taxonomy)

    DATA I/O OPTIONS:
        --outdir                <dir>   Final output dir for workflow results                   [default: 'results/nf_tram']
        --genome_fasta          <file>  Reference FASTA file for Trinity genome-guided assembly [default: unset]
        --entap_config_init     <file>  Input EnTap config file (only used as a starting point) [default: 'conf/entap_config.ini' in workflow repo]
        --entap_custom_db       <file>  FASTA file with a custom protein database               [default: none]
        --entap_nr_db           <file>  FASTA file with already downloaded NCBI NR database     [default: none => will be downloaded]
        --entap_refseq_db       <file>  FASTA file with already downloaded NCBI RefSeq database [default: none => will be downloaded]
        --entap_swissprot_db    <file>  FASTA file with already downloaded SwissProt database   [default: none => will be downloaded]
    
    MORE DATA I/O - PROVIDE WORKFLOW OUTPUTS TO SKIP STEPS:
        --genome_index          <dir>   STAR index dir for the '--genome_fasta' reference       [default: none => create index]
                                            This will cause the workflow to skip the reference indexing step.
        --genome_bam            <file>  BAM file of reads mapped to a ref. genome for genome-guided assembly with Trinity.
                                            Providing this file will make the workflow skip the map2genome steps.
        --reads_processed       <dir>   Path + glob pattern for with processed, non-normalized reads in FASTQ format.
        --reads_norm            <dir>   Directory with processed, normalized reads in FASTQ format.
        --assembly_1trans       <file>  Merged FASTA assembly with the longest transcript per gene (typically, EviGene output).
        --assembly_alltrans     <file>  Merged FASTA assembly with the all transcripts (typically, EviGene output).
        --entap_config_final    <file>  Final EnTap config file (after running EnTap config)
        --entap_diamond_dir     <dir>   Final EnTap DIAMOND db dir
        --entap_db_dir        #TODO EnTap binary DB etc, generated by config
        --assembly_final        <file>  Final FASTA assembly (typically, EnTap output).
        --tx2gene               <file>  A transcript-to-gene lookup file: a TSV with two columns. 

    OPTIONS TO DETERMINE WHAT PARTS OF THE WORKFLOW TO RUN:
        --subset_fq             <int>   Subset the FASTQ files to <int> reads                   [default: no subsetting]
        --skip_qc_reads         <bool>  Skip the read QC (FastQC => MultiQC) steps              [default: false]
        --skip_process_reads    <bool>  Skip the reads processing steps                         [default: false]
        --skip_assembly         <bool>  Skip the assembly steps                                 [default: false]
        --skip_assembly_norm
        --skip_assembly_nonorm
        --skip_spades
        --skip_trinity
        --skip_abyss
        --skip_qc_assembly      <bool>  Skip the assembly QC steps                              [default: false]
        --skip_annotate         <bool>  Skip the annotation steps                               [default: false]
        --skip_quantify         <bool>  Skip the expression quantification steps                [default: false]

    GENERAL SETTINGS:
        --strandedness          <str>   Sequencing library orientation: 'reverse', 'forward', or 'unstranded' [default: 'reverse']
        --trim_nextseq                  Use TrimGalore/Cutadapt's 'nextseq' option for poly-G trimming [default: don't use]
        --trim_qual             <int>   TrimGalore Phred min. base quality score                [default: 30]
        --trim_len              <int>   TrimGalore min. read length                             [default: 50]
        --nfiles_rcorr          <int>   Number of files to run Rcorrector with at a time        [default: 20 (=10 PE samples)]
        --kraken_db_url         <URL>   URL to a Kraken database/index from https://benlangmead.github.io/aws-indexes/k2
                                            [default: https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220926.tar.gz]
        --kraken_db_dir         <dir>   Path to a local Kraken database dir                     [default: none]
        --nfiles_assembly   #TODO
        --norms             #TODO
        --k_transabyss          <str>   Comma-separated list of kmer-values for TransAbyss      [default: '21,25,31,35,45,55,65,75,85']
        --k_spades              <str>   Comma-separated list of kmer-values for SPAdes          [default: '55,75,"auto"']
        --min_contig_length     <int>   Minimum contig length (for TransAbyss and Trinity)      [default: 300]
        --entap_qcov            <int>   EnTap min. query coverage for DIAMOND similarity search [default: 80]
        --entap_tcov            <int>   EnTap min. target coverage for DIAMOND similarity search[default: 80] 
        --entap_eval            <int>   EnTap evalue cutoff for DIAMOND similarity searching    [default: 1e-5]  
        --entap_fpkm            <num>   EnTap FPKM threshold for transcript filtering           [default: 0.5] 
        --entap_use_refseq #TODO
        --entap_use_nr          <bool>  Whether to use the NCBI NR database in EnTap            [default: true]
        --entap_use_tair        <bool>  Whether to use the (Arabidopsis) TAIR database in EnTap [default: true]         
        --entap_refseq_db_type  <str>   NCBI RefSeq database type for EnTap                     [default: 'complete']
                                            Options: 'complete', 'plant', 'vertebrate_mammalian', 'vertebrate_other', 'invertebrate'

    UTILITY OPTIONS
        --help                          Print this help message and exit
    """.stripIndent()
    exit 0
}


// =============================================================================
//                           GENERAL SETUP
// =============================================================================
// Print help if required
if (params.help) helpMessage()

// Include modules
include { SUBSET_FASTQ } from './modules/subset_fastq'

// Include subworkflows
include { QC_READS } from "./subworkflows/qc_reads" //addParams(OUTPUT: "${params.outdir}/read_qc")
include { PROCESS_READS } from "./subworkflows/process_reads"
include { NORM_READS } from "./subworkflows/norm_reads"
include { ASSEMBLE_NONORM } from "./subworkflows/assemble_nonorm"
include { ASSEMBLE_NORM } from "./subworkflows/assemble_norm"
include { MERGE_ASSEMBLIES } from "./subworkflows/merge_assemblies"
include { QC_ASSEMBLY } from "./subworkflows/qc_assembly"
include { ANNOTATE } from "./subworkflows/annotate"
include { QUANTIFY } from "./subworkflows/quantify"

// Check parameters
if (!params.reads) { exit 1, '\n============\nERROR: Input reads not specified! Use "--reads" to do so\n============\n' }
if (!params.busco_db) { exit 1, '\n============\nERROR: Busco db name not specified! Use "--busco_db" to do so\n============\n' }
if (!params.entap_taxon) { exit 1, '\n============\nERROR: Taxon name for EnTap not specified! Use "--entap_taxon" to do so\n============\n' }
if (!params.entap_contam) { exit 1, '\n============\nERROR: Contaminant list for EnTap not specified! Use "--entap_contam" to do so\n============\n' }

def checkPathParamList = [ params.reads, params.genome_index, params.genome_bam, params.reads_processed,
                            params.reads_norm, params.assembly_1trans, params.assembly_alltrans,
                            params.entap_config_final, params.entap_diamond_dir, params.entap_db_dir,
                            params.assembly_final, params.tx2gene ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Process paramaters
k_abyss = params.k_abyss?.split(',') as List // See https://github.com/nextflow-io/nextflow/discussions/2821
k_spades = params.k_spades?.split(',') as List
norms = params.norms?.split(',') as List


// =============================================================================
//                            DETERMINE WHAT TO RUN
// =============================================================================
// Process parameters
assemble_norm = params.skip_assembly_norm == true ? false : true
assemble_nonorm = params.skip_assembly_nonorm == true ? false : true
subset_fq = params.subset_fq == false ? false : true
subset_fq_count = params.subset_fq

// Starting points: run everything
qc_reads = true
process_reads = true
normalize = true
assemble = true
qc_assembly = true
annotate = true
quantify = true

// Determine what not to run, and load user-specified data
if (params.assembly_final != false && params.tx2gene != false) {
    qc_reads = false
    process_reads = false
    assemble = false
    qc_assembly = false
    annotate = false
    ch_assembly_final = Channel.fromPath(params.assembly_final, checkIfExists: true)
    ch_tx2gene = Channel.fromPath(params.tx2gene, checkIfExists: true)
}

if (params.assembly_1trans != false && params.assembly_alltrans != false && params.reads_processed != false) {
    qc_reads = false
    process_reads = false
    assemble = false
    ch_assembly_1trans = Channel.fromPath(params.assembly_1trans, checkIfExists: true)
    ch_assembly_alltrans = Channel.fromPath(params.assembly_alltrans, checkIfExists: true)
    ch_reads_processed = Channel.fromPath(params.reads_processed, checkIfExists: true)
}

if (params.reads_processed != false) {
    qc_reads = false
    process_reads = false
    ch_reads_processed = Channel.fromFilePairs(params.reads_processed, checkIfExists: true)
}

if (params.reads_norm != false) {
    normalize = false
    ch_reads_norm = Channel.fromFilePairs(params.reads_norm, checkIfExists: true)
    ch_reads_norm = ch_reads_norm.collect().flatten()
}

// Direct flags override run logic from above
if (params.skip_qc_reads == true) qc_reads = false
if (params.skip_process_reads == true) process_reads = false
if (params.skip_assembly == true) assemble = false
if (params.skip_qc_assembly == true) qc_assembly = false
if (params.skip_annotate == true) annotate = false
if (params.skip_quantify == true) quantify = false

// if 'assemble' is false, this overrides 'assemble_norm' and 'assemble_nonorm'
if (assemble == false) {
    assemble_norm = false
    assemble_nonorm = false
}


// =============================================================================
//                              REPORT
// =============================================================================
// Run info
log.info ""
log.info "======================================================================"
log.info "    T R A N S C R I P T O M E   A S S E M B L Y   P I P E L I N E"
log.info "======================================================================"
log.info "Parameters:"
params.each { k, v -> if (v) { println "${k}: ${v}" } }
log.info "======================================================================"
log.info ""

// Report what will be run
println "Subset FASTQ files?                    ${subset_fq}"
println "Run read QC?                           ${qc_reads}"
println "Run read processing?                   ${process_reads}"
println "Run transcriptome assembly?            ${assemble}"
println "   Run normalized assembly?            ${assemble_norm}"
println "   Run non-normalized assembly?        ${assemble_nonorm}"
println "Run assembly QC?                       ${qc_assembly}"
println "Run annotation?                        ${annotate}"
println "Run expression quantification?         ${quantify}"
log.info "======================================================================"


// =============================================================================
//                            THE WORKFLOW
// =============================================================================
// Default empty channels
ch_assemblies_norm = channel.empty()
ch_assemblies_nonorm = channel.empty()

// Workflow
workflow {
    ch_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
    println "Number of input samples:"; ch_reads.count().view()
    
    // Subset FASTQ files
    if (subset_fq) ch_reads = SUBSET_FASTQ(ch_reads, subset_fq_count).fq

    // Read QC
    if (qc_reads) QC_READS(ch_reads)
    
    // Read processing
    if (process_reads) {
        PROCESS_READS(ch_reads)
        ch_reads_processed = PROCESS_READS.out      // Format: Tuple by sample
    }

    if (assemble) {
        // Transcriptome assembly with normalized reads
        if (assemble_norm == true) {
            if (normalize == true) {
                NORM_READS(ch_reads_processed)
                ch_reads_norm = NORM_READS.out
            }
            ASSEMBLE_NORM(ch_reads_norm, k_abyss, k_spades)
            ch_assemblies_norm = ASSEMBLE_NORM.out
        }

        // Transcriptome assembly with non-normalized reads
        if (assemble_nonorm == true) {
            ASSEMBLE_NONORM(ch_reads_processed, k_abyss, k_spades)
            ch_assemblies_nonorm = ASSEMBLE_NONORM.out
        }

        // Merge assemblies
        ch_assemblies_all = ch_assemblies_norm.mix(ch_assemblies_nonorm).collect()
        MERGE_ASSEMBLIES(ch_assemblies_all)
        ch_assembly_1trans = MERGE_ASSEMBLIES.out.assembly_1trans
        ch_assembly_alltrans = MERGE_ASSEMBLIES.out.assembly_alltrans
    }
    
    // Assembly QC
    if (qc_assembly) {
        QC_ASSEMBLY(ch_assembly_1trans, ch_assembly_alltrans, ch_reads_processed)
    }

    // Assembly annotation
    if (annotate) {
        ANNOTATE(ch_assembly_1trans, ch_assembly_alltrans, ch_reads_processed)
        ch_assembly_final = ANNOTATE.out.assembly
        ch_tx2gene = ANNOTATE.out.tx2gene
    }

    // Transcript/gene quantification
    if (quantify) QUANTIFY(ch_assembly_final, ch_tx2gene, ch_reads_processed)
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
