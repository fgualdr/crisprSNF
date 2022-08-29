//
// Alignment with HISAT2
//

include { BWA_ALN       } from '../modules/bwa_align/align_main'
include { SAM_BAM_SORT_INDEX_SAMTOOLS       } from '../modules/samtools/sam_bam_sort_index'

workflow ALIGN_BWA {
    take:
    reads       // channel: [ val(meta), [ reads ] ]
    index       // channel: /path/to/star/index/

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with BWA_ALN
    //
    BWA_ALN ( reads, index  )
    ch_versions = ch_versions.mix(BWA_ALN.out.versions.first())

    // Substep 04.c
    // sam to BAM
    SAM_BAM_SORT_INDEX_SAMTOOLS (BWA_ALN.out.sam)
    ch_versions = ch_versions.mix(SAM_BAM_SORT_INDEX_SAMTOOLS.out.versions)

    emit:
    bam      = SAM_BAM_SORT_INDEX_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = SAM_BAM_SORT_INDEX_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats    = SAM_BAM_SORT_INDEX_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = SAM_BAM_SORT_INDEX_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAM_BAM_SORT_INDEX_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}