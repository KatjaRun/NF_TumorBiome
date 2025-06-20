/*
Processes to identify microorganisms in the fastq files containing the extracted unmapped reads
1) KRAKEN: Takes bam files with unmapped reads and assigns taxonomic labels, outputs .kreport2
2) BRACKEN: Takes .kreport2 and reeestimates the assigned counts based on genome size, outputs .bracken files
2) RECENTRIFUGE: 
3) BIOM: Collects all .bracken files and outputs them as a single biom file
*/

// Determine memory and set maxForks based on the calculated memory
include { calculateFolderSize; calculateMemoryForKraken } from '../conf/kraken_resources.groovy'
def krakenMemory = calculateMemoryForKraken(params.krakendb)
def krakenForks = (krakenMemory.toGiga() <= 100) ? 10 : 1 // Adjust maxForks based on memory (example logic)

process KRAKEN2 {
    tag "$barcode"
    memory { return krakenMemory } // Use pre-calculated memory
    maxForks krakenForks // Dynamically set maxForks
    
    publishDir "${params.outpath}/Kraken2"
    conda "./conf/Env.yml"


    input:
    tuple val(barcode), path(path_R1), path(path_R2)

    output:
    tuple val(barcode), path("${barcode}.kraken.txt"), path("${barcode}.kreport2")

    script:
    """
    # Running Kraken
    kraken2 \
        --db "${params.krakendb}" \
        --threads 20 \
        --output "./${barcode}.kraken.txt" \
        --paired "${path_R1}" "${path_R2}" \
        --report "./${barcode}.kreport2"
    """
}

process RECENTRIFUGE {
    tag "$platecenter"
    label "process_low", "error_ignore" // To ignore python write.xlsx error
    
    publishDir "${params.outpath}/Recentrifuge"
    conda "/home/rungger/.conda/envs/Recentrifuge"

    input:
    tuple val(platecenter), val(barcodes) 

    output:
    path "${platecenter}.rcf.data.tsv"
    path "${platecenter}.rcf.stat.tsv" 
    path "${platecenter}.rcf.html" 

    script:
    // Putting -k in front of every input to format it as input for Recentrifuge
    def kraken_format = barcodes.collect { "-k ${it}" }.join(' ')

    """
    rcf \
        ${kraken_format} \
        -n /home/rungger/.conda/envs/Recentrifuge/bin/taxdump \
        -o ${platecenter} \
        --scoring KRAKEN \
        --extra TSV \
        --minscore 35 \
        --exclude 9606 --takeoutroot 
    """
} 

process BRACKEN {
    tag "$barcode"
    label "standard"
    
    publishDir "${params.outpath}/Bracken"
    conda "./conf/Env.yml"


    input:
    tuple val(barcode), path(kraken_txt), path(kreport2)

    output:
    path("${barcode}.bracken")

    script:
    """
    # Running Bracken Reestimation
    bracken \
        -d "${params.krakendb}" \
        -i "${kreport2}" \
        -o "./${barcode}.bracken" \
        -l S
    """
}

process FINALBIOM {
    publishDir "${params.outpath}"
    conda "./conf/Env.yml"

    input:
    path bracken 

    output:
    path "Final.biom"

    script:
    """
    # Generating a final biom file for further analysis
    echo "Executing command: kraken-biom ${bracken} -o Final.biom" >&2
    kraken-biom ${bracken} -o Final.biom --fmt json
    """
}