//
// Alignment with HISAT2
//

include { BWAMEM_ALIGN       } from '../modules/bwa_align/align_main'
include { BAM_SORT_SAMTOOLS       } from '../modules/samtools/bam_sort_main'

workflow ALIGN_BWA {
    take:
    reads       // channel: [ val(meta), [ reads ] ]
    index       // channel: /path/to/star/index/

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with BWAMEM
    //
    BWAMEM_ALIGN ( reads, index )
    ch_versions = ch_versions.mix(BWAMEM_ALIGN.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_SAMTOOLS ( BWAMEM_ALIGN.out.bam )
    ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)

    emit:
    orig_bam = BWAMEM_ALIGN.out.bam           // channel: [ val(meta), bam   ]
    summary  = BWAMEM_ALIGN.out.summary       // channel: [ val(meta), log   ]
    fastq    = BWAMEM_ALIGN.out.fastq         // channel: [ val(meta), fastq ]

    bam      = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi      = BAM_SORT_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    stats    = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}