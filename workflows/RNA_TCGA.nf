#!/usr/bin/env nextflow

// Including all processes
include { GETBARCODE        } from '../modules/Preprocessing.nf'
include { AGGREGATEBARCODES } from '../modules/Preprocessing.nf'
include { FILTERUNMAPPED    } from '../modules/Preprocessing.nf'
include { KRAKEN2           } from '../modules/Identification.nf'
include { RECENTRIFUGE      } from '../modules/Identification.nf'
include { BRACKEN           } from '../modules/Identification.nf'
include { FINALBIOM         } from '../modules/Identification.nf'


// Workflow
workflow TCGARNA_identification {

    // Creating channel from either directory or samplesheet.csv
    if (params.input.endsWith('.csv')) {
        
        bam_ch = Channel.fromPath(params.input, checkIfExists: true)
            | splitCsv(header: true) 
            | map { row -> file(row.bam) }

    } else if (file(params.input).isDirectory()) {
        
        bam_ch = Channel.fromPath("${params.input}/*", type: "dir")
    
    } else {
        error "The input must either be a directoy containing BAM files or a CSV file."
    }
    
    // Extracting barcode and creating sample file
    sample_barcodes = GETBARCODE(bam_ch)
        | collect // Collecting all barcode files to get one Sample.tsv file
    final_barcodes = AGGREGATEBARCODES(sample_barcodes)

    // Creating a channel from the barcode output, with Barcode and path to BAM file
    barcode_ch = final_barcodes
        | splitCsv(sep:"\t", header: true) 
        | map { row -> [row.Barcode, file(row.Bam_Path)] }

    // Filtering out unmapped reads and getting tuple with barcode and R1 and R2 dirs
    unmapped_fastq = FILTERUNMAPPED(barcode_ch)

    // Kraken
    kraken = KRAKEN2(unmapped_fastq)

    kraken_ch = kraken
        | map { barcode, kraken_txt, kreport2 ->
        def last7 = barcode[-7..-1]
        return [last7, kraken_txt]
        }
        | groupTuple()

    // Running Recentrifuge
    RECENTRIFUGE(kraken_ch)

    // Bracken
    bracken = BRACKEN(kraken)
        .collect() //Collecting all bracken files to create one biom file

    // Create Biom file
    biom = FINALBIOM(bracken)

    // Decontam

    // Downstream

}