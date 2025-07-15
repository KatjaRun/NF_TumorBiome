#!/usr/bin/env nextflow

/*
========================================================================================
    NF_TumorBiome
========================================================================================
*/

nextflow.enable.dsl=2

// Printing out pipeline information
log.info """\
    NF_TumorBiome
    ============================================
    Metadata          : ${params.metadata}
    Abundance table   : ${params.abundance_table}
    Output            : ${params.outdir}
    ============================================
    """
    .stripIndent(true)

// Validating parameters with nf-schema
/*
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
validateParameters()
log.info paramsSummaryLog(workflow)
*/

// Including Workflows and Processes
include { CreateLineageTable         } from './workflows/Gathering_Data.nf'
include { Gathering_Data             } from './workflows/Gathering_Data.nf'
include { Basic_Analyses             } from './workflows/Basic_Analyses.nf'

/*
include { Transcript_counts_Analyses } from './workflows/Transcript_counts_Analyses.nf'
include { Survival_Analyses          } from './workflows/Survival_Analyses.nf'
*/

// Generating channels


workflow {

    // Creating channels of metadata, abundance table, and create lineage table if non-existent
    metadata_ch = channel.fromPath(params.metadata, checkIfExists: true)
    abundance_table_ch = channel.fromPath(params.abundance_table, checkIfExists: true)
    lineage_table_ch = params.lineage_table ? 
        channel.fromPath(params.lineage_table, checkIfExists: true) : 
        CreateLineageTable(abundance_table_ch)

    // Generating Phyloseq
    phyloseq_ch = Gathering_Data(metadata_ch, abundance_table_ch,lineage_table_ch)

    // Warnings if either samples or taxids were not completely matching
    phyloseq_ch.removed_samples.subscribe { path ->
        println "⚠️  Warning: Some samples were not matching. Removed samples are in: removed_samples.txt."
    }
    phyloseq_ch.removed_taxids.subscribe { path ->
        println "⚠️  Warning: Some taxids were not matching. Removed taxids are in: removed_taxids.txt."
    }

    // Basic Analyses
    Basic_Analyses(phyloseq_ch.phyloseq)
}

// Logging info for the end (hopefully)
workflow.onComplete {
    log.info "Pipeline completed successfully"
}
workflow.onError {
    log.error "Pipeline failed with an error"
}
