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
    path "Phyloseq.rds"

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
        // Build phyloseq
        physeq = CreatePhyloseq(metadata, abundance_table, lineage_table)

    emit:
        physeq
}