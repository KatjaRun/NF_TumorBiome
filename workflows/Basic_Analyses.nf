nextflow.enable.dsl=2

process Diversity {
    publishDir "${params.outdir}/Basic_Analyses/Diversity", mode: 'copy'

    input:
    path phyloseq
    
    output: 
    path "Alpha_Diversity.tsv"
    path "AlphaDiversity_*.png"
    path "BetaDiversity_NMDS.png"
    path "BetaDiversity_PCoA.png"

    script:
    """
    Diversity.R $phyloseq
    """
}


process Core {
    publishDir "${params.outdir}/Basic_Analyses/Core", mode: 'copy'

    input:
    path phyloseq
    
    output: 
    path "CoreParameters_*.tsv"
    path "CorePhyloseq_*.rds", optional: true, emit: core_rds
    path "CoreHeatmap_*.png", optional: true

    script:
    """
    Core.R $phyloseq
    """
}


process Spiec_Easi {
    label 'standard' 
    publishDir "${params.outdir}/Basic_Analyses/SPIEC_EASI", mode: 'copy'

    input:
    path core_phyloseq
    
    output: 
    path "Spiec_glasso_*.rds", optional: true

    script:
    """
    SpiecEasi.R $core_phyloseq
    """
}


workflow Basic_Analyses {

    take:
        phyloseq

    main:
        // Run diversity analyses
        Diversity(phyloseq)
        // Run core analyses and collect phyloseqs
        core_results = Core(phyloseq)
        core_phyloseqs = core_results.core_rds.flatten()
        // Run SpiecEasi on core phyloseqs
        Spiec_Easi(core_phyloseqs)
}
