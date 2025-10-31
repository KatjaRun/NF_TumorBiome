nextflow.enable.dsl=2

// To implement: also for unpaired data
process Kraken { 
    label 'standard'
    publishDir "${params.outdir}/Unmapped_Preprocessing/Kraken2", mode: 'copy'
    tag "$sample"

    input:
    tuple val(sample), path(unmapped_fastq_1), path(unmapped_fastq_2)

    output:
    path "${sample}.kraken.txt", emit: krakentxt
    path "${sample}.kreport2"

    script:
    """
    kraken2 \
        --db ${params.kraken_db} \
        --threads 20 \
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
    tuple val(sample), path(unmapped_fastq_1), path(unmapped_fastq_2)

    output:
    path "${sample}.results.txt", emit: centrifugetxt
    path "${sample}.report.txt"

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
    label 'standard'
    publishDir "${params.outdir}/Unmapped_Preprocessing/Recentrifuge", mode: 'copy'

    input:
    path classified
    val input_type
    val minscore

    output:
    path "All.rcf.data.tsv", emit: recentrifugetsv
    path "All.rcf.html"
    path "All.rcf.stat.tsv"

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
        -o All \
        --scoring ${scoring} \
        --extra TSV \
        --minscore ${minscore} \
        --exclude 9606 \
        --takeoutroot
    """
}

process Recentrifuge_to_abundance {
    label 'standard' 
    publishDir "${params.outdir}/Analysis_data", mode: 'copy'

    input:
    path rcf_file
    
    output: 
    path "Counts_rcf.tsv", emit: abundance_rcf

    script:
    """
    Rcf_to_abundance.R ${rcf_file}
    """
}


workflow Unmapped_Preprocessing {

    take:
    unmapped_ch

    main:
    if (params.use_centrifuge) {
        classified_centrifuge = Centrifuge(unmapped_ch)
        classified_centrifuge_ch = classified_centrifuge.centrifugetxt.collect().view()
        recentrifuge = Recentrifuge(classified_centrifuge_ch, 'centrifuge', params.rcf_minscore)
        recentrifuge_ch = recentrifuge.recentrifugetsv.flatten()
    } else {
        classified_kraken = Kraken(unmapped_ch) // default
        classified_kraken_ch = classified_kraken.krakentxt.collect().view()
        recentrifuge = Recentrifuge(classified_kraken_ch, 'kraken', params.rcf_minscore)
        recentrifuge_ch = recentrifuge.recentrifugetsv.flatten()
    }

    abundance_table = Recentrifuge_to_abundance(recentrifuge_ch)

    emit:
    abundance_table
}
