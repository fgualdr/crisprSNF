//
// Make index from fasta of gRNA using BWA
//

include { GUNZIP as GUNZIP_FASTA        } from '../modules/gunzip/gunzip'
include { UNTAR as UNTAR_BWA_INDEX      } from '../modules/untar/untar'
include { CUSTOM_GETCHROMSIZES          } from '../modules/getchromsizes/getchromsizes'
include { BWA_INDEX                     } from '../modules/bwa_index/bwa_index'

workflow PREPARE_GENOME {
    
    take:
    prepare_tool_indices // list   : tools to prepare indices for can be either BWA or Bowtie2

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = file(params.fasta)
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Uncompress BWA index or generate from scratch if required
    //
    ch_bwa_index = Channel.empty()
    if ('bwa' in prepare_tool_indices) {
        if (params.bwa_index) {
            if (params.bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( [ [:], params.bwa_index ] ).untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR_BWA_INDEX.out.versions)
            } else {
                ch_bwa_index = file(params.bwa_index)
            }
        } else {
            ch_bwa_index = BWA_INDEX ( ch_fasta ).index
            ch_versions     = ch_versions.mix(BWA_INDEX.out.versions)
        }
    }

    emit:
    fasta            = ch_fasta            //    path: genome.fasta
    fai              = ch_fai              //    path: genome.fai
    chrom_sizes      = ch_chrom_sizes      //    path: genome.sizes
    bwa_index        = ch_bwa_index     //    path: bwa/index/

    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}