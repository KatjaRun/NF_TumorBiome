nextflow.enable.dsl=2

// Checks if samples are matching and gets protein coding with gene name
process Check_Transcriptome {
    publishDir "${params.outdir}/HostTranscriptome_Analyses/", mode: 'copy'

    input:
    path phyloseq
    path host_transcriptome
    
    output: 
    path "HostTPM_filtered.tsv", emit: host_tpm
    path "Phyloseq_filtered.rds", emit: phyloseq_tpm
    path "Phyloseq_filtered_clr.rds", emit: phyloseq_tpm_clr
    path "removed_samples.txt", optional: true, emit: removed_samples

    script:
    """
    Check_HostTranscriptome.R $phyloseq $host_transcriptome
    """
}

// Running immune cell deconvolution and linear correlation
process Deconvolution {
    publishDir "${params.outdir}/HostTranscriptome_Analyses/Immune_Cells", mode: 'copy'

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

// Running ssGSEA and linear correlation
process Gsea {
    publishDir "${params.outdir}/HostTranscriptome_Analyses/Immune_Cells", mode: 'copy'

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

workflow HostTranscriptome_Analyses {

    take:
        phyloseq
        hosttranscriptome

    main:
        // Run diversity analyses
        Check_Transcriptome(phyloseq, hosttranscriptome)
        // Run core analyses and collect phyloseqs
        //Deconvolution(phyloseq)
        // Run SpiecEasi on core phyloseqs
        //Gsea(core_phyloseqs)

    emit: 
        removed_samples = Check_Transcriptome.out.removed_samples
}
