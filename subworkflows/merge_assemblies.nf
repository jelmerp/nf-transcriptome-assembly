// Subworkflow to merge multiple transcriptome assemblies

// Include modules
include { CONCAT_ASSEMBLIES } from '../modules/concat_assemblies'
include { EVIGENE } from '../modules/evigene'

// Define the subworkflow
workflow MERGE_ASSEMBLIES {
    take:
        Assemblies
    
    main:
        ch_assemblies_concat = CONCAT_ASSEMBLIES(Assemblies)
        ch_evigene = EVIGENE(ch_assemblies_concat.assembly)

    emit:
        assembly_1trans = ch_evigene.assembly_1trans
        assembly_alltrans = ch_evigene.assembly_alltrans
}
