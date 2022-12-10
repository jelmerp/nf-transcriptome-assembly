#!/usr/bin/env nextflow

// Nextflow transcriptome assembly pipeline

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
        --genome_index          <dir>   STAR index dir for the '--genome_fasta' reference       [default: none => create index]
                                            This will cause the workflow to skip the reference indexing step.
        --entap_config_init     <file>  Input EnTap config file (only used as a starting point) [default: 'conf/entap_config.ini' in workflow repo]
        --entap_custom_db       <file>  FASTA file with a custom protein database               [default: none]
        --entap_nr_db           <file>  FASTA file with already downloaded NCBI NR database     [default: none => will be downloaded]
        --entap_refseq_db       <file>  FASTA file with already downloaded NCBI RefSeq database [default: none => will be downloaded]
        --entap_swissprot_db    <file>  FASTA file with already downloaded SwissProt database   [default: none => will be downloaded]
    
    MORE DATA I/O - PROVIDE WORKFLOW OUTPUTS TO SKIP STEPS:
        --assembly_final        <file>  Final FASTA assembly (typically, EnTap output).
                                            Must be used with --tx2gene and --reads_nonorm. Only quantification will be run.
        --tx2gene               <file>  A transcript-to-gene map.
                                            Must be used with --assembly_final. Only quantification will be run.
        --assembly_1trans       <file>  Merged FASTA assembly with the longest transcript per gene (typically, EviGene output).
                                            Must be used with --assembly_alltrans and --reads_nonorm.
                                            Only annotation, assembly qc, and quantification will be run.
        --assembly_alltrans     <file>  Merged FASTA assembly with the all transcripts (typically, EviGene output).
                                            Must be used with --assembly_1trans and --reads_nonorm.
        --reads_norm            <dir>   Directory with processed, normalized reads in FASTQ format.
                                            Must be used with --reads_nonorm. Read QC and processing will be skipped.
        --reads_nonorm          <dir>   Directory with processed, non-normalized reads in FASTQ format.
        --entap_config_final    <file>  Final EnTap config file (after running EnTap config)
                                            Must be used with --entap_diamond_db_dir to skip the EnTap config step.
        --entap_diamond_db_dir  <dir>   Final EnTap DIAMOND db dir (after running EnTap config)
                                            Must be used with --entap_config_final to skip the EnTap config step.
        --genome_bam            <file>  BAM file of reads mapped to a ref. genome for genome-guided assembly with Trinity.
                                            Providing this file will make the workflow skip the map2genome steps.

    OPTIONS TO DETERMINE WHAT PARTS OF THE WORKFLOW TO RUN:
        --subset_fq             <int>   Subset the FASTQ files to <int> reads                   [default: no subsetting]
        --skip_qc_reads         <bool>  Skip the read QC (FastQC => MultiQC) steps              [default: false]
        --skip_process_reads    <bool>  Skip the reads processing steps                         [default: false]
        --skip_assembly         <bool>  Skip the assembly steps                                 [default: false]
        --skip_qc_assembly      <bool>  Skip the assembly QC steps                              [default: false]
        --skip_annotate         <bool>  Skip the annotation steps                               [default: false]
        --skip_quantify         <bool>  Skip the expression quantification steps                [default: false]

    GENERAL SETTINGS:
        --strandedness          <str>   Sequencing library orientation: 'reverse', 'forward', or 'unstranded' [default: 'reverse']
        --trim_nextseq                  Use TrimGalore/Cutadapt's 'nextseq' option for poly-G trimming [default: don't use]
        --trim_qual             <int>   TrimGalore Phred min. base quality score                [default: 5]
        --trim_len              <int>   TrimGalore min. read length                             [default: 36]
        --nfiles_rcorr          <int>   Number of files to run Rcorrector with at a time        [default: 20 (=10 PE samples)]
        --kraken_db_url         <URL>   URL to a Kraken database/index from https://benlangmead.github.io/aws-indexes/k2
                                            [default: https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220926.tar.gz]
        --kraken_db_dir         <dir>   Path to a local Kraken database dir                     [default: none]
        --k_transabyss          <str>   Comma-separated list of kmer-values for TransAbyss      [default: '21,25,31,35,45,55,65,75,85']
        --k_spades              <str>   Comma-separated list of kmer-values for SPAdes          [default: '55,75,"auto"']
        --min_contig_length     <int>   Minimum contig length (for TransAbyss and Trinity)      [default: 300]
        --entap_qcov            <int>   EnTap min. query coverage for DIAMOND similarity search [default: 80]
        --entap_tcov            <int>   EnTap min. target coverage for DIAMOND similarity search[default: 80] 
        --entap_eval            <int>   EnTap evalue cutoff for DIAMOND similarity searching    [default: 1e-5]  
        --entap_fpkm            <num>   EnTap FPKM threshold for transcript filtering           [default: 0.5] 
        --entap_use_nr          <bool>  Whether to use the NCBI NR database in EnTap            [default: true]
        --entap_use_tair        <bool>  Whether to use the (Arabidopsis) TAIR database in EnTap [default: true]         
        --entap_refseq_db_type  <str>   NCBI RefSeq database type for EnTap                     [default: 'complete']
                                            Options: 'complete', 'plant', 'vertebrate_mammalian', 'vertebrate_other', 'invertebrate'

    UTILITY OPTIONS
        --help                          Print this help message and exit
    """.stripIndent()
}

//  print help if required
if (params.help) {
    helpMessage()
    exit 0
}

// Check parameters
if (!params.reads) { exit 1, '\n============\nERROR: Input reads not specified! Use "--reads" to do so\n============\n' }
if (!params.busco_db) { exit 1, '\n============\nERROR: Busco db name not specified! Use "--busco_db" to do so\n============\n' }
if (!params.entap_taxon) { exit 1, '\n============\nERROR: Taxon name for EnTap not specified! Use "--entap_taxon" to do so\n============\n' }
if (!params.entap_contam) { exit 1, '\n============\nERROR: Contaminant list for EnTap not specified! Use "--entap_contam" to do so\n============\n' }

def checkPathParamList = [ params.reads ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Run info
log.info ""
log.info "======================================================================"
log.info "    T R A N S C R I P T O M E   A S S E M B L Y   P I P E L I N E"
log.info "======================================================================"
log.info "Parameters:"
params.each { k, v -> if (v) { println "${k}: ${v}" } }
log.info "======================================================================"
log.info ""

// Process paramaters
k_abyss = params.k_abyss?.split(',') as List // See https://github.com/nextflow-io/nextflow/discussions/2821
k_spades = params.k_spades?.split(',') as List

// Determine what to run
skip_qc_reads = false
skip_process_reads = false
skip_assembly = false
skip_qc_assembly = false
skip_annotate = false
skip_quantify = false

if (params.assembly_final != false && params.tx2gene != false) {
    skip_qc_reads = true
    skip_process_reads = true
    skip_assembly = true
    skip_qc_assembly = true
    skip_annotate = true
    assembly_final = Channel.fromPath(params.assembly_final)
    tx2gene = Channel.fromPath(params.tx2gene)
} else if (params.assembly_1trans != false && params.assembly_alltrans != false && params.reads_nonorm != false) {
    skip_qc_reads = true
    skip_process_reads = true
    skip_assembly = true
    assembly_1trans = Channel.fromPath(params.assembly_1trans)
    assembly_alltrans = Channel.fromPath(params.assembly_alltrans)
    reads_nonorm = Channel.fromPath(params.reads_nonorm)
} else if (params.reads_norm != false && params.reads_nonorm != false) {
    skip_qc_reads = true
    skip_process_reads = true
    reads_norm = Channel.fromPath(params.reads_norm)
    reads_nonorm = Channel.fromPath(params.reads_nonorm)
}

if (params.skip_qc_reads == true) skip_qc_reads = true
if (params.skip_process_reads == true) skip_process_reads = true
if (params.skip_assembly == true) skip_assembly = true
if (params.skip_qc_assembly == true) skip_qc_assembly = true
if (params.skip_annotate == true) skip_annotate = true
if (params.skip_quantify == true) skip_quantify = true

println "Skip read QC?                      ${skip_qc_reads}"
println "Skip read processing?              ${skip_process_reads}"
println "Skip read assembly?                ${skip_assembly}"
println "Skip annotation?                   ${skip_qc_assembly}"
println "Skip assembly QC?                  ${skip_annotate}"
println "Skip expression quantification?    ${skip_quantify}"

// Include modules
include { SUBSET_FASTQ } from './modules/subset_fastq'

// Include subworkflows
include { QC_READS } from "./subworkflows/qc_reads" //addParams(OUTPUT: "${params.outdir}/read_qc")
include { PROCESS_READS } from "./subworkflows/process_reads"
include { ASSEMBLY } from "./subworkflows/assembly"
include { QC_ASSEMBLY } from "./subworkflows/qc_assembly"
include { ANNOTATE } from "./subworkflows/annotate"
include { QUANTIFY } from "./subworkflows/quantify"

// Workflow
workflow {
    ch_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
    println "Number of input samples:"; ch_reads.count().view()
    
    // Subset FASTQ files
    if (params.subset_fq != false) ch_reads = SUBSET_FASTQ(ch_reads, params.subset_fq).fq

    // Read QC
    if (skip_qc_reads == false) QC_READS(ch_reads)
    
    // Read processing
    if (skip_process_reads == false) {
        PROCESS_READS(ch_reads)
        reads_norm = PROCESS_READS.out.reads_norm
        reads_nonorm = PROCESS_READS.out.reads_nonorm
    }

    // Transcriptome assembly
    if (skip_assembly == false) {
        ASSEMBLY(reads_norm, reads_nonorm, k_abyss, k_spades)
        assembly_1trans = ASSEMBLY.out.1trans
        assembly_alltrans = ASSEMBLY.out.alltrans
    }

    // Assembly QC
    if (skip_qc_assembly == false) {
        QC_ASSEMBLY(assembly_1trans, assembly_alltrans, reads_nonorm)
    }

    // Assembly annotation
    if (skip_annotate == false) {
        ANNOTATE(assembly_1trans, assembly_alltrans, reads_nonorm)
        assembly_final = ANNOTATE.out.assembly
        tx2gene = ANNOTATE.out.tx2gene
    }

    // Transcript/gene quantification
    if (skip_quantify == false) {
        QUANTIFY(assembly_final, tx2gene, reads_nonorm)
    }
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
