// Subworkflow to QC a transcriptome assembly

// Include modules
include { BUSCO } from '../modules/busco'
include { RNAQUAST } from '../modules/rnaquast'
include { DETONATE } from '../modules/detonate'

// Define the subworkflow
workflow QC_ASSEMBLY {
    take:
        Assembly_1trans
        Assembly_alltrans
        Reads_nonorm

    main:
        BUSCO(Assembly_1trans, params.busco_db)
        RNAQUAST(Assembly_alltrans)
        DETONATE(Assembly_alltrans, Reads_nonorm)
}
