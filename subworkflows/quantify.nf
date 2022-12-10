// Subworkflow to get gene expression counts from a transcriptome assembly

// Include modules
include { KALLISTO_INDEX } from '../modules/kallisto_index'
include { KALLISTO_QUANT } from '../modules/kallisto_quant'
include { TXIMPORT } from '../modules/tximport'

// Define the subworkflow
workflow QUANTIFY {
    take:
        Assembly
        Tx2gene
        Reads

    main:
        ch_kallisto_idx = KALLISTO_INDEX(Assembly)
        ch_kallisto = KALLISTO_QUANT(Reads, ch_kallisto_idx.index)
        TXIMPORT(ch_kallisto.counts_h5.collect(), Tx2gene)
}
