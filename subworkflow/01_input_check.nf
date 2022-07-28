//
// Check input samplesheet and get read channels
//

process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    01_check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}


// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id             = row.sample
    meta.grna_start     = row.grna_start
    meta.grnalength     = row.grnalength
    meta.r1_scaffold_trim    = row.r1_scaffold_trim
    meta.r1_scaffold_pos    = row.r1_scaffold_pos
    meta.umi_read            = row.umi_read
    meta.umi_pos             = row.umi_pos
    meta.umi_length          = row.umi_length
    meta.r2_scaffold_trim    = row.r2_scaffold_trim
    meta.r2_scaffold_pos    = row.r2_scaffold_pos


    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }

    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }
    
    fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    
    return fastq_meta
}
