nextflow.enable.dsl=2


process Kraken {
    label 'standard' 
    publishDir "${params.outdir}/Unmapped_Preprocessing/Kraken2", mode: 'copy'
    tag "$sample"

    input:
    tuple val(sample), path(unmapped_fastq_1), path(unmapped_fastq_2)

    output:
    path "${sample}.kraken.txt"
    path "${sample}.kreport2"

    script:
    """
    kraken2 \
        --db ${params.kraken_db} \
        --threads 20 \
        --paired ${unmapped_fastq_1} ${unmapped_fastq_2} \
        --report ${sample}.kreport2 \
        --output ${sample}.kraken.txt
    """
}

/*
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
*/

workflow Unmapped_Preprocessing {

    take:
    unmapped_ch

    main:
    if (params.use_centrifuge) {
        classified = Centrifuge(unmapped_ch)
        //refined = Recentrifuge_Centrifuge(classified)
    } else {
        classified = Kraken(unmapped_ch) // default
        //refined = Recentrifuge_Kraken(classified)
    }

    /*
    abundance_table = Recentrifuge_to_abundance(refined)

    emit:
    abundance_table
    */
}
