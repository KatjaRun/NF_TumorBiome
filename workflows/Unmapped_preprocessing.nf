nextflow.enable.dsl=2

// Determine memory and set maxForks based on the calculated memory for kraken
include { calculateFolderSize; calculateMemoryForKraken } from '../conf/kraken_resources.groovy'
def krakenMemory = calculateMemoryForKraken(params.kraken_db)
def krakenForks = (krakenMemory.toGiga() <= 100) ? 10 : 2 // Adjust maxForks based on memory

// To implement: also for unpaired data
process Kraken { 
    publishDir "${params.outdir}/Unmapped_Preprocessing/Kraken2", mode: 'copy'
    tag "$sample"
    memory { return krakenMemory } // Use pre-calculated memory
    maxForks krakenForks // Dynamically set maxForks

    input:
    tuple val(sample), path(unmapped_fastq_1), path(unmapped_fastq_2), val(batch)

    output:
    path "${sample}.kraken.txt"
    path "${sample}.kreport2"
    tuple val(sample), path("${sample}.kraken.txt"), val(batch), emit: krakentxt

    script:
    """
    kraken2 \
        --db ${params.kraken_db} \
        --threads $task.cpus \
        --paired ${unmapped_fastq_1} ${unmapped_fastq_2} \
        --report ${sample}.kreport2 \
        --output ${sample}.kraken.txt
    """
}

// To implement: also for unpaired data
process Centrifuge {
    label 'centrifuge_conda'
    publishDir "${params.outdir}/Unmapped_Preprocessing/Centrifuge", mode: 'copy'
    tag "$sample"

    input:
    tuple val(sample), path(unmapped_fastq_1), path(unmapped_fastq_2), val(batch)

    output:
    path "${sample}.results.txt"
    path "${sample}.report.txt"
    tuple val(sample), path("${sample}.results.txt"), val(batch), emit: centrifugetxt

    script:
    """
    mkdir ./temp
    export TMPDIR=\$PWD/tmp

    centrifuge \
        -x ${params.centrifuge_db} \
        --temp-directory ./temp \
        -p 8 \
        -1 ${unmapped_fastq_1} \
        -2 ${unmapped_fastq_2} \
        --report-file ${sample}.report.txt \
        -S ${sample}.results.txt
    """
}

// To implement: negative samples
process Recentrifuge {
    label 'long'
    publishDir "${params.outdir}/Unmapped_Preprocessing/Recentrifuge/$batch", mode: 'copy'
    tag "$batch"

    input:
    tuple val(batch), path(classified)
    val input_type
    val minscore

    output:
    path "${batch}.rcf.data.tsv", emit: recentrifugetsv
    path "${batch}.rcf.html"
    path "${batch}.rcf.stat.tsv"

    script:
    // get input files and scoring scheme
    def opts, scoring
    if (input_type == 'kraken') {
        opts    = classified.collect { "-k ${it}" }.join(' ')
        scoring = 'KRAKEN'
    } else if (input_type == 'centrifuge') {
        opts    = classified.collect { "-f ${it}" }.join(' ')
        scoring = 'SHEL'
    } else {
        throw new IllegalArgumentException("Unsupported input_type: ${input_type}")
    }

    """
    rcf \
        ${opts} \
        -n /home/rungger/.conda/envs/Recentrifuge/bin/taxdump \
        -o ${batch} \
        --scoring ${scoring} \
        --extra TSV \
        --minscore ${params.rcf_minscore} \
        --exclude 9606 \
        --takeoutroot
    """
}

process Recentrifuge_to_abundance {
    label 'standard' 
    publishDir "${params.outdir}/Analysis_data", mode: 'copy'

    input:
    path rcf_files
    
    output: 
    path "Counts_rcf.tsv", emit: abundance_rcf
    path "Rcftoabundance_RunInfo.txt"

    script:
    """
    Rcf_to_abundance.R ${rcf_files}
    """
}


workflow Unmapped_Preprocessing {

    take:
    unmapped_ch

    main:
    // Prepare channel and add batch if not available
    unmapped_norm_ch = unmapped_ch.map { row ->
        tuple(
            row.sample,
            file(row.unmapped_fastq_1),
            file(row.unmapped_fastq_2),
            row.batch ?: 'ALL'
        )
    }

    if (params.use_centrifuge) {
        classified_ch = Centrifuge(unmapped_norm_ch).centrifugetxt
        input_type = 'centrifuge'
    } else {
        classified_ch = Kraken(unmapped_norm_ch).krakentxt
        input_type = 'kraken'
    }

    // Grouping by batches for Recentrifuge
    classified_by_batch_ch = classified_ch
        .map { sample, file, batch -> tuple(batch, file) }
        .groupTuple()

    recentrifuge = Recentrifuge(
        classified_by_batch_ch,
        input_type,
        params.rcf_minscore
    )

    // Flatten Recentrifuge output
    recentrifuge_ch = recentrifuge.recentrifugetsv.collect().view()

    abundance_table = Recentrifuge_to_abundance(recentrifuge_ch)

    emit:
    abundance_table.abundance_rcf
}
