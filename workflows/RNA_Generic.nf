#!/usr/bin/env nextflow

// Including all processes
include { RUNNFCORE } from '../modules/Preprocessing.nf'

// Workflow
workflow GenericRNA_identification {

    // Running nf-core/rnaseq pipeline
    RUNNFCORE()

}