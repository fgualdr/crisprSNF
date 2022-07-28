process CUTADAP {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:4.1--py39hbf8eff0_0' :
        'quay.io/biocontainers/cutadapt:4.1--py39hbf8eff0_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trimmed.fastq.gz")      , emit: reads
    tuple val(meta), path("*.cutadapt.log")         , emit: log
    path "versions.yml"                              , emit: versions

    tuple val(meta), path("*.untrimmed.fastq.gz")   , emit: untrimmed

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (cores < 1) cores = 1
        if (cores > 4) cores = 4
    }

    // Added soft-links to original fastqs for consistent naming in MultiQC
    def prefix = task.ext.prefix ?: "${meta.id}"

    // define r1_adapter, r2_adapter
    def r1_scaffold = meta.r1_scaffold_trim
    def r2_scaffold = meta.r2_scaffold_trim

    def r1_scaffold_pos = meta.r1_scaffold_pos
    def r2_scaffold_pos = meta.r2_scaffold_pos

    def glength = meta.grnalength
    def umi_length = meta.umi_length

    if (r1_scaffold_pos == '3p') {
        if (r2_scaffold_pos == '3p') {
            println '3p R1 trimming and 3p R2 trimming'

            """
            [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
            [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz

            cutadapt \\
                    $args \\
                    --json=${prefix}.cutadapt.json \\
                    --cores=$cores \\
                    -m 1 \\
                    -a $r1_scaffold \\
                    -A $r2_scaffold \\
                    -o ${prefix}_1.trimmed.fastq.gz \\
                    -p ${prefix}_2.trimmed.fastq.gz \\
                    --untrimmed-output ${prefix}_1.untrimmed.fastq.gz \\
                    --untrimmed-paired-output ${prefix}_2.untrimmed.fastq.gz \\
                    ${prefix}_1.fastq.gz \\
                    ${prefix}_2.fastq.gz > ${prefix}.cutadapt.log
            
            
            cutadapt -l $glength  -o ${prefix}_1.trimmed_length.fastq.gz ${prefix}_1.trimmed.fastq.gz 
            rm ${prefix}_1.trimmed.fastq.gz 
            mv ${prefix}_1.trimmed_length.fastq.gz ${prefix}_1.trimmed.fastq.gz 

            cutadapt -l $umi_length  -o ${prefix}_2.trimmed_length.fastq.gz ${prefix}_2.trimmed.fastq.gz 
            rm ${prefix}_2.trimmed.fastq.gz 
            mv ${prefix}_2.trimmed_length.fastq.gz ${prefix}_2.trimmed.fastq.gz 
                    
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cutadapt: \$(echo \$(cutadapt --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            END_VERSIONS
            """
        
        } else if (r2_scaffold_pos == '5p') {
            println '3p R1 trimming and 5p R2 trimming'
            // 3' R1 trimming and 5' R2 trimming
        
            """
            [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
            [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz

            cutadapt \\
                    $args \\
                    --json=${prefix}.cutadapt.json \\
                    --cores=$cores \\
                    -m 1 \\
                    --length $glength \\
                    -a $r1_scaffold \\
                    -G $r2_scaffold \\
                    -o ${prefix}_1.trimmed.fastq.gz \\
                    -p ${prefix}_2.trimmed.fastq.gz \\
                    --untrimmed-output ${prefix}_1.untrimmed.fastq.gz \\
                    --untrimmed-paired-output ${prefix}_2.untrimmed.fastq.gz \\
                    ${prefix}_1.fastq.gz \\
                    ${prefix}_2.fastq.gz > ${prefix}.cutadapt.log

            cutadapt -l $glength  -o ${prefix}_1.trimmed_length.fastq.gz ${prefix}_1.trimmed.fastq.gz 
            rm ${prefix}_1.trimmed.fastq.gz 
            mv ${prefix}_1.trimmed_length.fastq.gz ${prefix}_1.trimmed.fastq.gz 

            cutadapt -l $umi_length  -o ${prefix}_2.trimmed_length.fastq.gz ${prefix}_2.trimmed.fastq.gz 
            rm ${prefix}_2.trimmed.fastq.gz 
            mv ${prefix}_2.trimmed_length.fastq.gz ${prefix}_2.trimmed.fastq.gz 

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cutadapt: \$(echo \$(cutadapt --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            END_VERSIONS
            """
        
        } 

    }else if (r1_scaffold_pos == '5p') {
        if (r2_scaffold_pos == '3p') {
            println '5p R1 trimming and 3p R2 trimming'
            // 3' R1 trimming and 3' R2 trimming

            """
            [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
            [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz

            cutadapt \\
                    $args \\
                    --json=${prefix}.cutadapt.json \\
                    --cores=$cores \\
                    -m 1 \\
                    --length $glength \\
                    -g $r1_scaffold \\
                    -A $r2_scaffold \\
                    -o ${prefix}_1.trimmed.fastq.gz \\
                    -p ${prefix}_2.trimmed.fastq.gz \\
                    --untrimmed-output ${prefix}_1.untrimmed.fastq.gz \\
                    --untrimmed-paired-output ${prefix}_2.untrimmed.fastq.gz \\
                    ${prefix}_1.fastq.gz \\
                    ${prefix}_2.fastq.gz > ${prefix}.cutadapt.log

            cutadapt -l $glength  -o ${prefix}_1.trimmed_length.fastq.gz ${prefix}_1.trimmed.fastq.gz 
            rm ${prefix}_1.trimmed.fastq.gz 
            mv ${prefix}_1.trimmed_length.fastq.gz ${prefix}_1.trimmed.fastq.gz 

            cutadapt -l $umi_length  -o ${prefix}_2.trimmed_length.fastq.gz ${prefix}_2.trimmed.fastq.gz 
            rm ${prefix}_2.trimmed.fastq.gz 
            mv ${prefix}_2.trimmed_length.fastq.gz ${prefix}_2.trimmed.fastq.gz 

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cutadapt: \$(echo \$(cutadapt --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            END_VERSIONS
            """

        
        } else if (r2_scaffold_pos == '5p') {
            println '5p R1 trimming and 5p R2 trimming'
            // 3' R1 trimming and 5' R2 trimming

            """
            [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
            [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz

            cutadapt \\
                    $args \\
                    --json=${prefix}.cutadapt.json \\
                    --cores=$cores \\
                    -m 1 \\
                    --length $glength \\
                    -g $r1_scaffold \\
                    -G $r2_scaffold \\
                    -o ${prefix}_1.trimmed.fastq.gz \\
                    -p ${prefix}_2.trimmed.fastq.gz \\
                    --untrimmed-output ${prefix}_1.untrimmed.fastq.gz \\
                    --untrimmed-paired-output ${prefix}_2.untrimmed.fastq.gz \\
                    ${prefix}_1.fastq.gz \\
                    ${prefix}_2.fastq.gz > ${prefix}.cutadapt.log

            cutadapt -l $glength  -o ${prefix}_1.trimmed_length.fastq.gz ${prefix}_1.trimmed.fastq.gz 
            rm ${prefix}_1.trimmed.fastq.gz 
            mv ${prefix}_1.trimmed_length.fastq.gz ${prefix}_1.trimmed.fastq.gz 

            cutadapt -l $umi_length  -o ${prefix}_2.trimmed_length.fastq.gz ${prefix}_2.trimmed.fastq.gz 
            rm ${prefix}_2.trimmed.fastq.gz 
            mv ${prefix}_2.trimmed_length.fastq.gz ${prefix}_2.trimmed.fastq.gz 

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cutadapt: \$(echo \$(cutadapt --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            END_VERSIONS
            """
        
        } 
    }




    
}
