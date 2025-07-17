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
    path "CoreHeatmap_*.png", optional: true

    script:
    """
    Core.R $phyloseq
    """
}

/*
process Spiec_Easi {

}
*/

workflow Basic_Analyses {

    take:
        phyloseq

    main:
        Diversity(phyloseq)
        Core(phyloseq)

}
