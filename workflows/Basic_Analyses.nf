nextflow.enable.dsl=2

process Diversity {
    label 'standard'
    publishDir "${params.outdir}/Basic_Analyses/Diversity", mode: 'copy'

    input:
    path phyloseq
    
    output: 
    path "AlphaDiversity.tsv", emit: phyloseq
    path "AlphaDiversity_*.png"
    path "BetaDiversity_NMDS.tsv"
    path "BetaDiversity_NMDS_*.png"
    path "BetaDiversity_PCoA.tsv"
    path "BetaDiversity_PCoA_*.png"
    path "Composition_*.png"
    path "Diversity_RunInfo.txt"

    script:
    """
    Diversity.R $phyloseq
    """
}


process Core {
    label 'standard'
    publishDir "${params.outdir}/Basic_Analyses/Core", mode: 'copy'

    input:
    path phyloseq
    
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
