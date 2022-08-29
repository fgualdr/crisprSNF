//
// Sam to BAM, Sort, index BAM file and run samtools stats, flagstat and idxstats
//

include { SAMTOOLS_SAM_BAM      } from '../../modules/samtools/sam_bam/samtools_sam_bam'
include { SAMTOOLS_SORT      } from '../../modules/samtools/sort/samtools_sort'
include { SAMTOOLS_INDEX     } from '../../modules/samtools/index/samtools_index'
include { BAM_STATS_SAMTOOLS } from './bam_stats_samtools'

workflow SAM_BAM_SORT_INDEX_SAMTOOLS {
    take:
    ch_sam // channel: [ val(meta), [ sam ] ]

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_SAM_BAM ( ch_sam )
    ch_versions = ch_versions.mix(SAMTOOLS_SAM_BAM.out.versions.first())

    SAMTOOLS_SORT ( SAMTOOLS_SAM_BAM.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .map {
            meta, bam, bai  ->
                if (bai) {
                    [ meta, bam, bai ]
                } 
        }
        .set { ch_bam_bai }

    BAM_STATS_SAMTOOLS ( ch_bam_bai )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]

    stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}