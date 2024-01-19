process MAGECK_TEST {
    tag "$meta.name"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mageck:0.5.9--py37h6bb024c_0':
        'biocontainers/mageck:0.5.9--py37h6bb024c_0' }"

    input:
    tuple val(meta), path(count_table)

    output:
    path("*.gene_summary.txt")  , emit: gene_summary
    path("*.sgrna_summary.txt") , emit: sgrna_summary
    path("*.R")                 , emit: r_script
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'mageck_test'
    def target = "${meta.target}"
    def control = "${meta.control}"
    def name = "${meta.name}"

    """
    mageck  \\
        test \\
        $args \\
        -k $count_table \\
        -t $target \\
        -c $control \\
        -n $name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mageck: \$(mageck -v)
    END_VERSIONS
    """
}
