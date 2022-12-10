// Subworkflow to annotate a transcriptome assembly (merged by EviGene)

// Parameters assumed to be passed from main.nf:
// - entap_config_init
// ... plus see below

// Process params
refseq_type = params.entap_refseq_type
refseq_db = params.entap_refseq_db
nr_db = params.entap_nr_db
use_nr_db = params.entap_use_nr_db
use_tair_db = params.entap_use_tair_db
swissprot_db = params.entap_swissprot_db
custom_db = params.entap_custom_db

if (params.entap_config_final != false && params.entap_diamond_db_dir != false) {
    ch_db_dir = Channel.fromPath(params.entap_diamond_db_dir)
    ch_final_config = Channel.fromPath(params.entap_config_final)
    skip_entap_config = true
    skip_get_protein_dbs = true
}

// Include modules
include { GET_REFSEQ } from '../modules/get_refseq'
include { GET_NR } from '../modules/get_nr'
include { GET_TAIR } from '../modules/get_tair'
include { MAP2TRANSCRIPTOME } from '../modules/map2transcriptome'
include { ENTAP_CONFIG } from '../modules/entap_config'
include { ENTAP } from '../modules/entap'
include { ENTAP_PROCESS } from '../modules/entap_process'

// Define the subworkflow
workflow ANNOTATE {
    take:
        Assembly_1trans     // Transcriptome assembly with only the longest transcript per gene
        Assembly_alltrans   // Transcriptome assembly with all transcripts
        Reads               // Non-normalized, processed reads in FASTQ format

    main:
        // Map reads to the transcriptome
        ch_map2trans = MAP2TRANSCRIPTOME(Assembly_1trans, Reads)

        // Get protein reference databases
        if (skip_get_protein_dbs == false) {
            ch_refseq = refseq_db == false ? GET_REFSEQ(refseq_type) : channel.fromPath(refseq_db)
            ch_tair = use_tair_db == true ? GET_TAIR() : channel.empty()
            if (use_nr_db == true) {
                ch_nr = nr_db == false ? GET_NR() : channel.fromPath(nr_db)
            } else {
                ch_nr = channel.empty()
            }
            ch_custom_db = custom_db ? channel.fromPath(refseq_db) : channel.empty()
            ch_ref_dbs_fa = ch_refseq.mix(ch_nr, ch_tair, ch_custom_db).collect()
        }

        // Configure and run EnTap
        if (skip_entap_config == false) {
            ch_entap_conf = ENTAP_CONFIG(params.entap_config_init, ch_ref_dbs_fa)
            ch_db_dir = ch_entap_conf.db_dir
            ch_final_config = ch_entap_conf.config
        }
        
        ch_entap = ENTAP(Assembly_1trans, ch_db_dir, ch_final_config, ch_map2trans.bam)
        
        ch_entap_ed = ENTAP_PROCESS(ch_entap.out, Assembly_1trans, Assembly_alltrans)

    emit:
        assembly = ch_entap_ed.assembly
        tx2gene = ch_entap_ed.tx2gene
}
