nextflow.enable.dsl=2

process Run_nfcore_rnaseq {
    input:
    path samplesheet

    output:
    path "results"

    """
    nextflow run nf-core/rnaseq \
        --input ${samplesheet} \
        --outdir results \
        --save_unaligned \
        --skip_alignment \
        -profile docker
    """
}

process Kraken {

}

process Centrifuge {

}

process Recentrifuge {
    input:
    path classified_report

    output:
    path "recentrifuge_output/*" emit: out

    script:
    def report_type = params.use_centrifuge ? "centrifuge" : "kraken2"
    
    """
    mkdir -p recentrifuge_output
    recentrifuge -i ${classified_report} -o recentrifuge_output --type ${report_type}
    """
}


process Recentrifuge_to_abundance {

}

workflow Rnaseq_Preprocessing {

    take:
    fastq

    main:
    if (params.use_nfcore) {
        preprocessed = Run_nfcore_rnaseq(samplesheet: file(params.samplesheet_rnaseq))
        unaligned = preprocessed.out
    } else {
        unaligned = fastq
    }

    if (params.use_centrifuge) {
        classified = Centrifuge(unaligned)
        refined = Recentrifuge_Centrifuge(classified)
    } else {
        classified = Kraken(unaligned) // default
        refined = Recentrifuge_Kraken(classified)
    }

    abundance_table = Recentrifuge_to_abundance(refined)

    emit:
    abundance_table
}
