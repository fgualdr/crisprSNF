process BWA_ALN {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7' :
        'quay.io/biocontainers/bwa:0.7.17--hed695b0_7' }"

    input:
    tuple val(meta), path(reads)
    path index


    output:
    tuple val(meta), path("*.sam"), emit: sam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    bwa aln \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        > ${prefix}.sai
        
    bwa samse \\
        \$INDEX \\
        ${prefix}.sai \\
        $reads \\
        > ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

}

