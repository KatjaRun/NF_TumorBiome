/*
Processes to preprocess TCGA Bam files to prepare for Kraken and Bracken
1) GETBARCODE: Gets TCGA Barcode from the file UUID and creates a Barcode.tsv file
2) AGGREGATEBARCODES: Collects all Barcode.tsv files from previous process and produces one Sample.tsv file
3) FILTERUNMAPPED: filters out unmapped reads from BAM files and creates two FASTQ files with quality filtered reads for R1 and R2
4) GATHERBATCHES
*/


process GETBARCODE {
    maxForks 20 // Limiting number of tasks to avoid accession error
    
    conda "./conf/Env.yml"
    
    input: 
    val bamdir 

    output: 
    path "*.tsv"

    script:
    """ 
    Rscript ${baseDir}/bin/Get_Barcode.R '${bamdir}'
    """
}


process AGGREGATEBARCODES {
    publishDir params.outpath, mode: 'copy' 

    input:
    path sample_barcodes

    output:
    path "Sample_File.tsv"

    script:
    def cat_files = sample_barcodes.collect{ "cat $it >> Sample_File.tsv" }.join("\n") // Function to collect all files and write them to Samples.tsv
    """
    # Giving the file a header
    echo -e "Barcode\tPlateCenter\tUUID\tBam_Path" >> Sample_File.tsv
    # Writing the individual lines to the file
    ${cat_files}
    """
}

process FILTERUNMAPPED {
    tag "$barcode"
    label "process_low"
    maxForks 10
    
    publishDir "${params.outpath}/Fastqs"
    conda "./conf/Env.yml"
    
    input:
    tuple val(barcode), path(input_bam)  // Expect a tuple with barcode and BAM file path

    output:
    tuple val(barcode), path("${barcode}_R1.fq.gz"), path("${barcode}_R2.fq.gz")

    script:
    """
    # Run the sh script to extract unmapped reads
    Filter_Unmapped.sh '${barcode}' '${input_bam}'
    """
}

process RUNNFCORE {
    label "error_ignore"
    
    script:
    """
    nextflow run nf-core/rnaseq -r 3.14.0 \
    --input ${params.input} \
    --outdir ${params.outpath} \
    --kraken_db ${params.krakendb} \
    --save_kraken_assignments true \
    --save_kraken_unassigned true \
    --gtf /data/genomes/hg38/annotation/gencode/gencode.v44.primary_assembly.annotation.gtf \
    --fasta /data/genomes/hg38/fasta/gencode/GRCh38.primary_assembly.genome.fa \
    --star_index "/data/genomes/hg38/index/STAR/2.7.9a/gencode/gencode.v44.GRCh38.primary_assembly.genome/200" \
    --skip_bigwig true \
    --skip_trimming false \
    --gencode true \
    --aligner "star_salmon" \
    -profile singularity,icbi \
    -resume
    """
}

/*
For 3.17.0
nextflow run nf-core/rnaseq \
--input /data/projects/2020/OvarianCancerHH/Data_FFPE_Pancreas/Nextflow_Results/Samplesheet.csv \
--outdir /data/projects/2020/OvarianCancerHH/Data_FFPE_Pancreas/Nextflow_Results \
--kraken_db /data/databases/kraken2/k2_pluspf_20240605 \
--save_kraken_assignments true \
--save_kraken_unassigned true \
--gtf /data/genomes/hg38/annotation/gencode/gencode.v44.primary_assembly.annotation.gtf \
--fasta /data/genomes/hg38/fasta/gencode/GRCh38.primary_assembly.genome.fa \
--star_index "/data/genomes/hg38/index/STAR/2.7.9a/gencode/gencode.v44.GRCh38.primary_assembly.genome/200" \
--skip_bigwig true \
--skip_trimming false \
--gencode true \
--aligner "star_salmon" \
-profile singularity,icbi \
-resume
-work-dir

For 3.14.0
nextflow run nf-core/rnaseq -r 3.14.0 \
--input /data/projects/2020/OvarianCancerHH/Data_FFPE_Pancreas/Nextflow_Results/Samplesheet.csv \
--outdir /data/projects/2020/OvarianCancerHH/Data_FFPE_Pancreas/Nextflow_Results \
--gtf /data/genomes/hg38/annotation/gencode/gencode.v44.primary_assembly.annotation.gtf \
--fasta /data/genomes/hg38/fasta/gencode/GRCh38.primary_assembly.genome.fa  \
--star_index /data/genomes/hg38/index/STAR/2.7.9a/gencode/gencode.v44.GRCh38.primary_assembly.genome/200 \
--skip_bigwig true \
--skip_trimming false \
--gencode true \
--aligner star_salmon \
--save_unaligned true \
-profile singularity,icbi \
-resume \
-work-dir /data/projects/2020/OvarianCancerHH/Data_FFPE_Pancreas/Nextflow_Results/work

*/

