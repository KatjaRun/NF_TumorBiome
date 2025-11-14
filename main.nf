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
    Samplesheet       : ${params.samplesheet}
    Metadata          : ${params.metadata}
    Abundance table   : ${params.abundance_table}
    Host transcriptome: ${params.host_transcriptome}
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
include { Unmapped_Preprocessing     } from './workflows/Unmapped_preprocessing.nf'
include { CreateLineageTable         } from './workflows/Gathering_Data.nf'
include { Gathering_Data             } from './workflows/Gathering_Data.nf'
include { Basic_Analyses             } from './workflows/Basic_Analyses.nf'
include { HostTranscriptome_Analyses } from './workflows/HostTranscriptome_Analyses.nf'

/*
include { Survival_Analyses          } from './workflows/Survival_Analyses.nf'
*/

// Generating channels


workflow {

    // Generating channel for metadata (is needed as sample names need to be provided)
    metadata_ch = Channel.fromPath(params.metadata, checkIfExists: true)

    if (params.samplesheet) {
        // Generate abundance and lineage table if unmapped preprocessing is performed
        unmapped_ch = Channel.fromPath(params.samplesheet)
                            .splitCsv(header: true)

        abundance_table_ch = Unmapped_Preprocessing(unmapped_ch)
        lineage_ch = CreateLineageTable(abundance_table_ch)

    } else if (params.abundance_table) {
        abundance_table_ch = Channel.fromPath(params.abundance_table, checkIfExists: true)
        lineage_ch = params.lineage_table ?
            Channel.fromPath(params.lineage_table, checkIfExists: true) :
            CreateLineageTable(abundance_table_ch)
    } else {
        error "You must provide either a samplesheet (--samplesheet) or an abundance table (--abundance_table)."
    }
    
    // Generating Phyloseq
    phyloseq_ch = Gathering_Data(metadata_ch, abundance_table_ch,lineage_ch)

    // Warnings if either samples or taxids were not completely matching
    phyloseq_ch.removed_samples.subscribe { path ->
        println "⚠️  Warning: Some samples were not matching. Removed samples are in: Analysis_data/removed_samples.txt."
    }
    phyloseq_ch.removed_taxids.subscribe { path ->
        println "⚠️  Warning: Some taxids were not matching. Removed taxids are in: Analysis_data/removed_taxids.txt."
    }

    // Basic Analyses
    Basic_Analyses(phyloseq_ch.phyloseq)

    // Host-transcriptome analyses if host_transcriptome is provided
    if (params.host_transcriptome) {
        ht_ch = Channel.fromPath(params.host_transcriptome, checkIfExists: true)
        ht_checked_ch = HostTranscriptome_Analyses(phyloseq_ch.phyloseq, ht_ch)

        // Warning if samples between tpm and phyloseq were not completely matching
        ht_checked_ch.removed_samples.subscribe { path ->
            println "⚠️  Warning: Some samples were not matching. Removed samples are in: HostTranscriptome_Analyses/removed_samples.txt."
        }
    }
}


// Logging info for the end (hopefully)
workflow.onComplete {
    log.info "Pipeline completed successfully"
}
workflow.onError {
    log.error "Pipeline failed with an error"
}
