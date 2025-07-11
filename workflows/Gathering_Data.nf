nextflow.enable.dsl=2

process CreateLineageTable {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path abundance_table
    
    output: 
    path "Lineage_table.csv"
    

    script:
    """
    Taxonomy.py $abundance_table
    """
}

process CreatePhyloseq {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path metadata
    path abundance_table
    path lineage_table

    output:
    path "Phyloseq.rds", emit: phyloseq
    path "removed_samples.txt", optional: true
    path "removed_taxids.txt", optional: true

    script:
    """
    Create_Phyloseq.R $metadata $abundance_table $lineage_table
    """
}

workflow Gathering_Data {

    take:
        metadata
        abundance_table
        lineage_table

    main:
        CreatePhyloseq(metadata, abundance_table, lineage_table)

    emit:
        CreatePhyloseq.out.phyloseq
}
