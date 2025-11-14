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
    path phyloseq_clr
    path host_transcriptome
    
    output: 
    path "Cibersort.tsv"
    path "Quantiseq.tsv"
    path "Microbiome_clr_*.tsv"
    path "Spearman_Cibersort_*.tsv"
    path "Spearman_Quantiseq_*.tsv"
    path "Heatmap_Cibersort_*.png"
    path "Heatmap_Quantiseq_*.png"

    script:
    """
    Deconvolution.R $phyloseq_clr $host_transcriptome
    """
}

// Running ssGSEA and linear correlation
process Gsea {
    publishDir "${params.outdir}/HostTranscriptome_Analyses/Immune_Parameters", mode: 'copy'

    input:
    path phyloseq_clr
    path host_transcriptome
    
    output: 
    path "GSEA_Reactome.tsv"
    path "GSEA_GOBP.tsv"
    path "GSEA_Hallmarks.tsv"
    path "Spearman_Reactome_*.tsv"
    path "Spearman_GOBP_*.tsv"
    path "Spearman_Hallmarks_*.tsv"
    path "Heatmap_Reactome_*.png"
    path "Heatmap_GOBP_*.png"
    path "Heatmap_Hallmarks.png"

    script:
    """
    Gsea.R $phyloseq_clr $host_transcriptome
    """
}

workflow HostTranscriptome_Analyses {

    take:
        phyloseq
        hosttranscriptome

    main:
        // Run diversity analyses
        checkt_results = Check_Transcriptome(phyloseq, hosttranscriptome)
        // Run core analyses and collect phyloseqs
        Deconvolution(checkt_results.phyloseq_tpm_clr, checkt_results.host_tpm)
        // Run SpiecEasi on core phyloseqs
        Gsea(checkt_results.phyloseq_tpm_clr, checkt_results.host_tpm)

    emit: 
        removed_samples = Check_Transcriptome.out.removed_samples
}
