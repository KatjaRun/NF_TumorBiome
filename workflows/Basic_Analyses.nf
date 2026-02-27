nextflow.enable.dsl=2

process Diversity {
    label 'standard'
    tag "$name"
    publishDir "${params.outdir}/Basic_Analyses/Diversity/$name", mode: 'copy'

    input:
    tuple val(name), path(phyloseq)
    
    output: 
    path "AlphaDiversity_*.tsv", emit: phyloseq
    path "AlphaDiversity_*.png"
    path "BetaDiversity_NMDS_*.tsv", optional: true
    path "BetaDiversity_NMDS_*.png", optional: true
    path "BetaDiversity_PCoA.tsv", optional: true
    path "BetaDiversity_PCoA_*.png", optional: true
    path "Composition_*.png"
    path "Diversity_RunInfo.txt"

    script:
    """
    Diversity.R $phyloseq $name
    """
}


process Core {
    label 'standard'
    tag "$name"
    publishDir "${params.outdir}/Basic_Analyses/Core/$name", mode: 'copy'

    input:
    tuple val(name), path(phyloseq)
    
    output: 
    path "CoreParameters_*.tsv"
    path "CorePhyloseq_*.rds", optional: true, emit: core_rds
    path "CoreHeatmap_*.png", optional: true
    path "Core_RunInfo.txt"

    script:
    """
    Core.R $phyloseq
    """
}


process Spiec_Easi {
    label 'standard' 
    publishDir "${params.outdir}/Basic_Analyses/SPIEC_EASI", mode: 'copy'
    tag "$core_phyloseq"

    input:
    path core_phyloseq
    
    output: 
    path "Spiec_glasso_*.rds", optional: true
    //path "SpiecEasi_RunInfo_*.txt"

    script:
    """
    SpiecEasi.R $core_phyloseq
    """
}


workflow Basic_Analyses {

    take:
        phyloseq_all

    main:
        // Run diversity analyses
        Diversity(phyloseq_all)
        // Run core analyses and collect phyloseqs
        core_results = Core(phyloseq_all)
        core_phyloseqs = core_results.core_rds.flatten()
        // Run SpiecEasi on core phyloseqs
        Spiec_Easi(core_phyloseqs)
}
