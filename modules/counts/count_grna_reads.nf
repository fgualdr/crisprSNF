process COUNT_GRNA_READS {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3' :
        'quay.io/biocontainers/bedtools:2.30.0--h468198e_3' }"

    input:
    path(bams)
    path(bai)
    path(bed)

    output:
    path("*counts_per_gRNA.txt")        , emit: counts
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    println "bams $bams "
    println "bed $bed "

    """
    bedtools multicov \\
        $args \\
        -bams ${bams.join(' ')} \\
        -bed $bed > tmp.txt

    awk '
    BEGIN { FS = OFS = "\\t" ; print "chr", "start", "end", "${bams.join('","')}" }
    { print \$0, "" }
    ' tmp.txt > counts_per_gRNA.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1)
    END_VERSIONS
    """
}