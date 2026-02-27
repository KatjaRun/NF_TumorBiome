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
    path "Phyloseq_Bacteria.rds", optional: true, emit: phyloseq_bacteria
    path "Phyloseq_Viruses.rds", optional: true, emit: phyloseq_virus
    path "Phyloseq_Archaea.rds", optional: true, emit: phyloseq_archaea
    path "Phyloseq_Fungi.rds", optional: true, emit: phyloseq_fungi
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
        cp = CreatePhyloseq(metadata, abundance_table, lineage_table)

    emit:
        phyloseq_mic    = cp.phyloseq_microbes
        phyloseq_all = Channel.empty()
            .concat(cp.phyloseq_microbes.map { file -> tuple('Microbes', file) })
            .concat(cp.phyloseq_bacteria.map { file -> tuple('Bacteria', file) })
            .concat(cp.phyloseq_virus.map { file -> tuple('Viruses', file) })
            .concat(cp.phyloseq_archaea.map { file -> tuple('Archaea', file) })
            .concat(cp.phyloseq_fungi.map { file -> tuple('Fungi', file) })
        removed_samples = cp.removed_samples
        removed_taxids  = cp.removed_taxids
}
