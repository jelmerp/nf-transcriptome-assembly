// Subworkflow to annotate a transcriptome assembly (merged by EviGene)

// Include modules
include { GET_REFSEQ } from '../modules/get_refseq'
include { GET_NR } from '../modules/get_nr'
include { GET_TAIR } from '../modules/get_tair'
include { DIAMOND_MAKEDB } from '../modules/diamond_makedb'
include { MAP2TRANSCRIPTOME } from '../modules/map2transcriptome'
include { ENTAP_CONFIG } from '../modules/entap_config'
include { ENTAP } from '../modules/entap'
include { ENTAP_PROCESS } from '../modules/entap_process'
include { FOFN } from '../modules/fofn'

// Pre-existing files
nr_db = params.entap_nr_db
refseq_db = params.entap_refseq_db
custom_db = params.entap_custom_db
refseq_type = params.entap_refseq_type
swissprot_db = params.entap_swissprot_db
config_init = params.entap_config_init
config_final = params.entap_config_final
diamond_dir = params.entap_diamond_dir

// Determine what to run
filter_expression = params.entap_expr_filt
use_refseq = params.entap_use_refseq
use_nr = params.entap_use_nr
use_tair = params.entap_use_tair
run_entap_config = true
get_protein_dbs = true

if (config_final != false && diamond_dir != false) {
    ch_diamond_dir = Channel.fromPath(diamond_dir)
    ch_config_final = Channel.fromPath(config_final)
    run_entap_config = false
    get_protein_dbs = false
}

// Default empty channels
ch_nr = channel.empty()
ch_refseq = channel.empty()
ch_tair = channel.empty()
ch_custom_db = channel.empty()
ch_bam = channel.empty()

ch_config_init = channel.fromPath(config_init)

// Define the subworkflow
workflow ANNOTATE {
    take:
        Assembly_1trans     // Transcriptome assembly with only the longest transcript per gene
        Assembly_alltrans   // Transcriptome assembly with all transcripts
        Reads               // Non-normalized, processed reads in FASTQ format

    main:
        // Map reads to the transcriptome
        if (filter_expression) {
            ch_read_fofn = FOFN(Reads.collect().flatten().collect(), "fastq.gz")
            ch_bam = MAP2TRANSCRIPTOME(Assembly_1trans, ch_read_fofn).bam
        }

        // Get protein reference databases
        if (get_protein_dbs) {
            // RefSeq
            if (use_refseq) ch_refseq = refseq_db == false ? GET_REFSEQ(refseq_type) : channel.fromPath(refseq_db)
            
            // NR
            if (use_nr) ch_nr = nr_db == false ? GET_NR() : channel.fromPath(nr_db)
            
            // TAIR (Arabidopsis)
            if (use_tair) ch_tair = GET_TAIR() 

            // Custom DB
            if (custom_db != false) ch_custom_db = channel.fromPath(custom_db)
            
            // Combine all FASTAs into a channel
            ch_ref_dbs_fa = ch_tair.mix(ch_nr, ch_refseq, ch_custom_db).collect()
            
            // Create DIAMOND dbs from protein FASTAs
            ch_diamond = DIAMOND_MAKEDB(ch_ref_dbs_fa)
            ch_diamond = ch_diamond.db.collect()
        }

        // Configure and run EnTap
        if (run_entap_config) ch_config = ENTAP_CONFIG(ch_config_init).config
        ch_entap = ENTAP(Assembly_1trans, ch_config, ch_diamond, ch_bam).out
        ch_entap_ed = ENTAP_PROCESS(ch_entap, Assembly_1trans, Assembly_alltrans)

    emit:
        assembly = ch_entap_ed.assembly
        tx2gene = ch_entap_ed.tx2gene
}
