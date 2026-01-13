nextflow.enable.dsl=2

process CreateLineageTable {
    label 'standard' 
    publishDir "${params.outdir}/Analysis_data", mode: 'copy'

    input:
    path abundance_table
    
    output: 
    path "Lineage_table.tsv"
    

    script:
    """
    Taxonomy.R $abundance_table
    """
}

process CreatePhyloseq {
    label 'standard' 
    publishDir "${params.outdir}/Analysis_data", mode: 'copy'

    input:
    path metadata
    path abundance_table
    path lineage_table

    output:
    path "Phyloseq_Full.rds", emit: phyloseq
    path "Phyloseq_Microbes.rds", emit: phyloseq_microbes
    path "Phyloseq_*.rds", emit: phyloseq_specific
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
        phyloseq        = CreatePhyloseq.out.phyloseq_microbes
        removed_samples = CreatePhyloseq.out.removed_samples
        removed_taxids  = CreatePhyloseq.out.removed_taxids
}
