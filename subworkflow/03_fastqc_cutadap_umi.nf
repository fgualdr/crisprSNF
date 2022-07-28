//
// Read QC, RM scaffold, select gRNA and UMI 
//

include { FASTQC           } from '../modules/fastqc/main'
include { CUTADAP       } from '../modules/cutadap/main'
include { UMIEXTRACT       } from '../modules/umitools/main'

// include { UMITOOLS_EXTRACT } from '../modules/umitools/main'

workflow FASTQC_CUTADAP_UMI {

    take:
    reads            // channel: [ val(meta), [ reads ] ]

    main:

    // i) fastqc before the trimming
    ch_versions = Channel.empty()
    fastqc_html = Channel.empty()
    fastqc_zip  = Channel.empty()

    FASTQC ( reads ).html.set { fastqc_html }
    fastqc_zip  = FASTQC.out.zip
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    // ii) Cupadapers 
    trim_reads    = reads
    untrimmed = Channel.empty()
    trim_log      = Channel.empty()
    
    CUTADAP ( reads ).reads.set { trim_reads }
    untrimmed = CUTADAP.out.untrimmed
    trim_log = CUTADAP.out.log
    ch_versions   = ch_versions.mix(CUTADAP.out.versions.first())

    // iii) assign Umi to R1 discard R2
    umi_reads    = trim_reads
    umi_log      = Channel.empty()

    UMIEXTRACT ( trim_reads ).reads.set { umi_reads }
    umi_log = CUTADAP.out.log

    emit:

    reads = umi_reads // channel: [ val(meta), [ reads ] ]

    untrimmed  = untrimmed    // channel: [ val(meta), [ reads ] ]
    trim_log   = trim_log        // channel: [ val(meta), [ log ] ]
    umi_log   = umi_log        // channel: [ val(meta), [ log ] ]

    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
    
}
