nextflow.enable.dsl=2

process CreateLineageTable {
    publishDir "${params.outdir}/Analysis_data", mode: 'copy'

    input:
    path abundance_table
    
    output: 
    path "Lineage_table.tsv"
    

    script:
    """
    Taxonomy.py $abundance_table
    """
}

process CreatePhyloseq {
    publishDir "${params.outdir}/Analysis_data", mode: 'copy'

    input:
    path metadata
    path abundance_table
    path lineage_table

    output:
    path "Phyloseq.rds", emit: phyloseq
    path "removed_samples.txt", optional: true, emit: removed_samples
    path "removed_taxids.txt", optional: true, emit: removed_taxids

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
        phyloseq        = CreatePhyloseq.out.phyloseq
        removed_samples = CreatePhyloseq.out.removed_samples
        removed_taxids  = CreatePhyloseq.out.removed_taxids
}
