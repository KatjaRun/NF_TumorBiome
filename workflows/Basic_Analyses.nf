nextflow.enable.dsl=2

process Diversity {
    publishDir "${params.outdir}", mode: 'copy'

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

/*
process Core {

}

process Spiec_Easi {

}
*/

workflow Basic_Analyses {

    take:
        phyloseq

    main:
        Diversity(phyloseq)

}
