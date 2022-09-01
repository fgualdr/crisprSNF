//
// Read QC, RM scaffold, select gRNA and UMI 
//

include { FASTQC           } from '../modules/fastqc/fastqc'
include { CUTADAP       } from '../modules/cutadap/cutadap'
include { UMIEXTRACT       } from '../modules/umitools/umiextract'

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

    // iii) assign Umi to R1 discard R2
    umi_reads    = reads
    umi_log      = Channel.empty()

    UMIEXTRACT ( reads ).reads.set { umi_reads }
    umi_log = UMIEXTRACT.out.log

    // ii) Cupadapers 
    trim_reads    = umi_reads
    untrimmed = Channel.empty()
    trim_log      = Channel.empty()
    
    CUTADAP ( umi_reads ).reads.set { trim_reads }
    untrimmed = CUTADAP.out.untrimmed
    trim_log = CUTADAP.out.log
    ch_versions   = ch_versions.mix(CUTADAP.out.versions.first())

    emit:

    reads = trim_reads // channel: [ val(meta), [ reads ] ]

    untrimmed  = untrimmed    // channel: [ val(meta), [ reads ] ]
    trim_log   = trim_log        // channel: [ val(meta), [ log ] ]
    umi_log   = umi_log        // channel: [ val(meta), [ log ] ]
    fastqc_zip = fastqc_zip
    
    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
    
}
