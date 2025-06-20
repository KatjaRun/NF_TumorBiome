#!/usr/bin/env nextflow

/*
NF Pipeline to identify microorganisms from bam files
*/

nextflow.enable.dsl=2

/* 
Parameter values
*/
params.input = null
params.outpath = null
params.krakendb = "/data/databases/kraken2/k2_pluspf_20240605"
params.type = null

/*
Printing out pipeline information
*/
log.info """\
    MICROBIAL IDENTIFICATION NF-PIPELINE
    ====================================
    Input       : ${params.input}
    Output path : ${params.outpath}
    Kraken DB   : ${params.krakendb}
    Running for : ${params.type}
    ====================================
    """
    .stripIndent(true)

/* 
Starting Workflow
*/
include { TCGARNA_identification }    from './workflows/RNA_TCGA.nf'
include { GenericRNA_identification } from './workflows/RNA_Generic.nf'

workflow {
        if (params.type == "genericRNA") {
        GenericRNA_identification()
    } else if (params.type == "TCGARNA") {
        TCGARNA_identification()
    } else {
        log.error "Invalid params.type value: '${params.type}'. It should be 'genericRNA' or 'TCGRNA'."
        exit 1
    }
}
