#!/usr/bin/env nextflow

nextflow.enable.dsl = 2 

/*
========================================================================================
                         crisprSNF
========================================================================================
 Analysis Pipeline.
 #### https://github.com/fgualdr/crisprSNF
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:
      nextflow run ./crisprSNF --input design.csv --fasta ... -profile singularity
    Mandatory arguments:
      --input [file]                  Comma-separated file containing information about the samples in the experiment (see docs/usage.md) (Default: './design.csv')
      --fasta [file]                  Path to Fasta reference. This is the gRNA library in fasta format.
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, docker, singularity, awsbatch, test
    Other
      --outdir [file]                 The output directory where the results will be saved (Default: './results')
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic (Default: false)
    """.stripIndent()
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                SET UP CONFIGURATION VARIABLES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input,
    params.fasta,
    params.bwa_index,
    params.bowtie_index
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, 'Samples design file not specified!' }
if (params.aligner) { ch_aligner = params.aligner } else { exit 1, 'Aligner not specified! Can be bowtie2 or bwa' }
if (params.sort_bam) { ch_sort_bam = params.sort_bam } else { ch_sort_bam = true }
if (params.umitools_dedup_stats) { ch_umitools_dedup_stats = params.umitools_dedup_stats } else { ch_umitools_dedup_stats = true }

// Multi QC
ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
// Header files for MultiQC
ch_spp_nsc_header           = file("$projectDir/assets/multiqc/spp_nsc_header.txt", checkIfExists: true)
ch_spp_rsc_header           = file("$projectDir/assets/multiqc/spp_rsc_header.txt", checkIfExists: true)
ch_spp_correlation_header   = file("$projectDir/assets/multiqc/spp_correlation_header.txt", checkIfExists: true)
ch_peak_count_header        = file("$projectDir/assets/multiqc/peak_count_header.txt", checkIfExists: true)
ch_frip_score_header        = file("$projectDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
ch_peak_annotation_header   = file("$projectDir/assets/multiqc/peak_annotation_header.txt", checkIfExists: true)
ch_deseq2_pca_header        = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_deseq2_clustering_header = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)

def prepareToolIndices  = []
prepareToolIndices << params.aligner

if (params.fasta) {
    lastPath = params.fasta.lastIndexOf(File.separator)
    bwa_base = params.fasta.substring(lastPath+1)
    ch_fasta = file(params.fasta, checkIfExists: true)
} else {
    exit 1, 'Fasta file not specified!'
}

if (params.bwa_index) {
    lastPath = params.bwa_index.lastIndexOf(File.separator)
    bwa_dir  = params.bwa_index.substring(0,lastPath+1)
    bwa_base = params.bwa_index.substring(lastPath+1)
    Channel
        .fromPath(bwa_dir, checkIfExists: true)
        .set { ch_bwa_index }
}

if (params.bowtie_index) {
    lastPath = params.bowtie_index.lastIndexOf(File.separator)
    bwa_dir  = params.bowtie_index.substring(0,lastPath+1)
    bwa_base = params.bowtie_index.substring(lastPath+1)
    Channel
        .fromPath(bwa_dir, checkIfExists: true)
        .set { ch_bowtie_index }
}

def multiqc_report      = []

if (params.rm_samp) { ch_rm_samp = params.rm_samp } else {ch_rm_samp = 'NULL' }

println "input $params.input "
println "outdir $params.outdir "
println "publish_dir_mode $params.publish_dir_mode "
println "aligner $params.aligner "

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     LOAD SUBWORKFLOWS                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

include { INPUT_CHECK } from './subworkflow/01_input_check' 
include { CONCAT_FASTQ } from './subworkflow/02_cat_fastq'
include { FASTQC_CUTADAP_UMI } from './subworkflow/03_fastqc_cutadap_umi'
include { PREPARE_GENOME } from './subworkflow/04_prepare_genoma'
include { ALIGN_BWA } from './subworkflow/05_align_bwa'
include { DEDUP_UMI_UMITOOLS } from './subworkflow/06_dedup_umi_tools'
include { SAMTOOLS_SORT               } from './modules/samtools/sort/samtools_sort'
include { COUNT_GRNA_READS               } from './modules/counts/count_grna_reads'
include { NORMDEG } from './subworkflow/07_norm_deg'
include { MAGECK_TEST               } from './modules/mageck/test/'
include { MAGECK_MLE               } from './modules/mageck/mle/'
include { MULTIQC } from './modules/multiqc'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from './modules/custom/dumpsoftwareversions/main'

workflow CRISPRSNF {

    ch_versions = Channel.empty()

    //
    // STEP01: Read in samplesheet, validate and stage input files
    // This step takes as input a comma separated csv file including rid, sid, sample name id, single_paired, the path to the run to use, the input is required and the lanes to use
    // lanes has being included given that occasionally issues associated to sequencer lane might need to be expluded
    //

    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            [ meta_clone, fastq ] 
    }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1             // here it refers to single fastq per sample
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1              // here it refers to multiple fastq per sample
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    ch_versions.view { "value: $it" } 

    //
    // STEP 02: Concatenate FastQ files from same sample if required 
    // Based on the previous step concatenate fastqs per sample
    // 

    CONCAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CONCAT_FASTQ.out.versions.first().ifEmpty(null))

    //
    // STEP 03: 
    // Cutadapt 
    // Identify and remove adapter sequence from both R1 and R2.
    // Report identified position
    //

    ch_fastqc_raw_multiqc  = Channel.empty()
    ch_fastqc_trim_multiqc = Channel.empty()
    ch_trim_log_multiqc    = Channel.empty()

    FASTQC_CUTADAP_UMI (
        ch_cat_fastq
    )
    ch_versions = ch_versions.mix(FASTQC_CUTADAP_UMI.out.versions)
    ch_filtered_reads = FASTQC_CUTADAP_UMI.out.reads
    
    //
    // STEP 04: 
    // Mapping
    // This step include the creation of the index and the mapping
    // we test BWA-MEM and in future Bowtie2 in SE mode acccording to param.aligner

    // Check type of aligner

    if( ch_aligner == 'bwa' ){
        println "aligner bwa"

        if (!params.bwa_index){
            // STEP 04.a
            // Check BWA index
            PREPARE_GENOME (
                prepareToolIndices
            )
            ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

        }

        // STEP 04.b
        // Check align
        ALIGN_BWA (
            ch_filtered_reads,
            PREPARE_GENOME.out.bwa_index
        )

        ch_genome_bam        = ALIGN_BWA.out.bam
        ch_genome_bam_index  = ALIGN_BWA.out.bai
        ch_samtools_stats    = ALIGN_BWA.out.stats
        ch_samtools_flagstat = ALIGN_BWA.out.flagstat
        ch_samtools_idxstats = ALIGN_BWA.out.idxstats

        ch_versions = ch_versions.mix(ALIGN_BWA.out.versions)


    }else{
        exit 1, "ERROR: aligner must be bwa"
    }

    // STEP 05
    // DEDUP UMITOOLS
    DEDUP_UMI_UMITOOLS (
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                params.umitools_dedup_stats
    )
    ch_genome_bam        = DEDUP_UMI_UMITOOLS.out.bam
    ch_genome_bam_index  = DEDUP_UMI_UMITOOLS.out.bai

    ch_dedup_samtools_stats    = DEDUP_UMI_UMITOOLS.out.stats
    ch_dedup_samtools_flagstat = DEDUP_UMI_UMITOOLS.out.flagstat
    ch_dedup_samtools_idxstats = DEDUP_UMI_UMITOOLS.out.idxstats

    ch_versions = ch_versions.mix(DEDUP_UMI_UMITOOLS.out.versions)

    
    // STEP 06
    // Assemble a teble with gRNA x Samples :
    ch_genome_bam
            .collect{it[1]}
            .set { ch_counts_bam }
    ch_genome_bam_index
            .collect{it[1]}
            .set { ch_counts_bai }

    COUNT_GRNA_READS (
        ch_counts_bam,
        ch_counts_bai,
        PREPARE_GENOME.out.chrom_bed
    )
    ch_versions = ch_versions.mix(COUNT_GRNA_READS.out.versions)


    // STEP 07 Normalize and perform DEGs calling # R code in bin
    // we do this by defining the conditions which are the meta.id without the replicate i.e. the ending _R[1-9]


    ch_genome_bam
        .map { meta, bam ->
            def new_meta = [:]
            new_meta.id = meta.id
            new_meta.condition = meta.id.replaceAll(/^(.*)_R[1-9].*$/, '$1')
            [new_meta.condition]
        }
        .unique()
        .set { ch_conditions }

    ch_conditions
        .combine(ch_conditions)
        .filter { pair ->
            pair[0] != pair[1] && pair[0] != 'PLASMID' && pair[1] != 'PLASMID'
        }
        .map { pair ->
            def new_meta = [:]
            new_meta.name = "${pair[0]}_vs_${pair[1]}"
            new_meta.target = pair[0]
            new_meta.control = pair[1]
            [new_meta]
        }
        .set { ch_combined_meta_deg }

    ch_combined_meta_deg
        .combine(COUNT_GRNA_READS.out.counts)
        .set { ch_deg_test }
        
    ch_deg_test.view()

    // STEP 07 Normalize and perform DEGs calling # R code in bin
    // Perform only if DEG_DESIGN is present
    NORMDEG(
        ch_deg_test
    )

    // // STEP 09 MAGEK TEST
    // // GENERATE a tuple of comparisons
    // input is a meta - with count data .. 
    // We define paired conditions i.e. the count table has to be added to a tuple with -t and -c
    // We go by replicate so each column "PRESORT_R1" will be paired with any other column like "HIGH_R1" or "LOW_R1"
    // while "PRESORT_R2" will be paired with any other column like "HIGH_R2" or "LOW_R2" etc..
    // so we group samples by replicate
    // then we pair each sample with any other sample in the same replicate. generating a new meta.control and meta.treatment
    // generate a new channel where meta is modified.
    // we want the channel to include only the meta.id and meta.replicate and meta.condition
    // we get the meta.id from the ch_genome_bam
    // meta.replicate is generated by removing any character before the last '_R' in the meta.id
    // meta.condition is generated by removing any character after the last '_R' in the meta.id including the _R itself
    ch_genome_bam
        .map { meta, bam ->
            def new_meta = [:]
            new_meta.id = meta.id
            new_meta.replicate = (meta.id =~ /_R(\d+)/) ? "R${meta.id.replaceAll(/.*_R(\d+).*/, '$1')}" : "R0"
            new_meta.condition = meta.id.replaceAll(/^(.*)_R[1-9].*$/, '$1') 
            [new_meta]
        }
        .set { ch_meta }

    // This returns something like: [[id:RETNLA_LOW_R1, replicate:_R1, condition:LOW]]
    // Group elements by replicate and define meta.target and meta.control then combine with COUNT_GRNA_READS.out.counts_mageck
    // the end results should be [[meta], counts]
    
    ch_meta
        .combine(ch_meta)
        .filter { pair ->
            pair[0].replicate == pair[1].replicate && pair[0].condition != pair[1].condition
        }
        .map { pair ->
            def new_meta = [:]
            new_meta.name = "${pair[0].id}_vs_${pair[1].id}"
            new_meta.target = pair[0].id
            new_meta.control = pair[1].id
            [new_meta]
        }
        .set { ch_combined_meta }

    // ch_combined_meta is a new meta with the form: [[target:RETNLA_LOW_R1, control:RETNLA_HIGH_R1]]
    // each need to be paired with the unique counts table - so we need to repeat the counts table for each meta to get [[target:RETNLA_LOW_R1, control:RETNLA_HIGH_R1], counts]

    ch_combined_meta
        .combine(COUNT_GRNA_READS.out.counts_mageck)
        .set { ch_mageck_test }


    // we group by replicate and define meta.target and meta.control
    // the global design of the assay include samples labelled PRESORT_R[0-9], HIGH_R[0-9], LOW_R[0-9], PLASMID_R[0-9]
    // basically we want to generate achannel of meta.target and meta.control where each meta.target is paired with each meta.control by replicate.
    // Targets are HIGH_R[0-9], LOW_R[0-9] and the PRESORT_R[0-9] is the control
    // we group by replicate and define meta.target and meta.control


    MAGECK_TEST(
         ch_mageck_test
     )


    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

    // MULTIQC (
    //         ch_multiqc_config,
    //         ch_multiqc_custom_config.collect().ifEmpty([]),
    //         CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),

    //         ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_trim_log_multiqc.collect{it[1]}.ifEmpty([]),

    //         ch_samtools_stats.collect{it[1]}.ifEmpty([]),
    //         ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
    //         ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),

    //         ch_dedup_samtools_stats.collect{it[1]}.ifEmpty([]),
    //         ch_dedup_samtools_flagstat.collect{it[1]}.ifEmpty([]),
    //         ch_dedup_samtools_idxstats.collect{it[1]}.ifEmpty([]),

    //         MARK_DUPLICATES_PICARD.out.stats.collect{it[1]}.ifEmpty([]),
    //         MARK_DUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]),
    //         MARK_DUPLICATES_PICARD.out.idxstats.collect{it[1]}.ifEmpty([]),
    //         MARK_DUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]),

    //         BAM_FILTER_EM.out.stats.collect{it[1]}.ifEmpty([]),
    //         BAM_FILTER_EM.out.flagstat.collect{it[1]}.ifEmpty([]),
    //         BAM_FILTER_EM.out.idxstats.collect{it[1]}.ifEmpty([]),
    //         ch_picardcollectmultiplemetrics_multiqc.collect{it[1]}.ifEmpty([]),
    
    //         ch_deeptoolsplotprofile_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_deeptoolsplotfingerprint_multiqc.collect{it[1]}.ifEmpty([]),
    
    //         ch_custompeaks_frip_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_custompeaks_count_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_plothomerannotatepeaks_multiqc.collect().ifEmpty([]),
    //         ch_subreadfeaturecounts_multiqc.collect{it[1]}.ifEmpty([]),

    //         ch_deseq2_pca_multiqc.collect().ifEmpty([]),
    //         ch_deseq2_clustering_multiqc.collect().ifEmpty([])
    //     )
    //     multiqc_report = MULTIQC.out.report.toList()

}


workflow {
    CRISPRSNF ()
}