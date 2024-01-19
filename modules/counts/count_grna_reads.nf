process COUNT_GRNA_READS {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.0--hf5e1c6e_3' :
        'quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_3 ' }"

    input:
    path(bams)
    path(bai)
    path(bed)

    output:
    path("counts_per_gRNA.txt")            , emit: counts
    path("counts_for_gRNA_mageck.txt")     , emit: counts_mageck
    path "versions.yml"                     , emit: versions

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

    awk '\\
    BEGIN { FS = OFS = "\\t" ; print "chr", "start", "end", "${bams.join('","')}" }\\
    { print \$0, "" }' tmp.txt > tmp_counts_per_gRNA.txt

    # replace from the first line of tmp_counts_per_gRNA.txt the ".bam" by nothing and save to file counts_for_gRNA_mageck
    sed '1s/.bam//g' tmp_counts_per_gRNA.txt > counts_per_gRNA.txt
    rm tmp_counts_per_gRNA.txt

    # we modify the tmp.txt file to be compatible with mageck
    awk -F'\\t' 'BEGIN {OFS = "\\t"} {\\
        split(\$1, a, /_/);\\
        for (i = 1; i <= 2; i++) {\\
            printf \$1 "\\t" a[1] "\t";\\
            for (j = 4; j <= NF; j++) {\\
                printf \$j;\\
                if (j < NF) {\\
                    printf "\\t";\\
                } else {\\
                    printf "\\n";\\
                }\\
            }\\
        }\\
    }' tmp.txt > tmp_mageck.txt

    awk '\\
    BEGIN { FS = OFS = "\\t" ; print "sgRNA", "Gene", "${bams.join('","')}" }\\
    { print \$0, "" }\\
    ' tmp_mageck.txt > counts_for_gRNA_mageck_a.txt

    # replace from the first line of counts_for_gRNA_mageck_a.txt the ".bam" by nothing and save to file counts_for_gRNA_mageck
    sed '1s/.bam//g' counts_for_gRNA_mageck_a.txt > counts_for_gRNA_mageck.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1)
    END_VERSIONS
    """
}