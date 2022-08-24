process UMIEXTRACT {
    tag "$meta.id"
    label "process_low"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.2--py38hbff2b2d_1' :
        'quay.io/biocontainers/umi_tools:1.1.2--py38hbff2b2d_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def umi_length = meta.umi_length.toInteger()
    def pattern = 'N'*umi_length

    // For now the library uses UMI only on one R2
    // For this design we consider this the only
    // as we use umi_tools it expect the UMI always on the primary read
    // Therefore we need to swap r1 and r2:

    """
    umi_tools extract \\
            -I ${reads[1]} \\
            --read2-in=${reads[0]} \\
            -S ${prefix}.umi_extract_2.fastq.gz \\
            --read2-out=${prefix}.umi_extract_1.fastq.gz \\
            --bc-pattern=$pattern \\
            $args \\
            > ${prefix}.umi_extract.log

    rm -rf ${prefix}.umi_extract_2.fastq.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            umitools: \$(umi_tools --version 2>&1 | sed 's/^.*UMI-tools version://; s/ *\$//')
    END_VERSIONS
    """
    
}