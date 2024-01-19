process NORMDEG {
    tag "$meta.name"
    label 'process_high_memory'
    container 'docker://fgualdr/envnorm'
    echo true

    input:
    tuple val(meta), path(count_table)

    output:
    path('*CrisprRes*')     , emit: path_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def skip_outliers = params.skip_outliers ? params.skip_outliers : false
    def target = "${meta.target}"
    def control = "${meta.control}"
    def name = "${meta.name}"
    
    """

    02_norm_diff_kernelTarget.r \\
            -f ${count_table} \\
            -p ${task.cpus} \\
            -o ${skip_outliers} \\
            -t ${target} \\
            -c ${control}
    """
    
}
