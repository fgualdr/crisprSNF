#!/usr/bin/env nextflow

nextflow.enable.dsl = 2 

/*
========================================================================================
                         crisprSNF
========================================================================================
 Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/
----------------------------------------------------------------------------------------
*/

def helpMessage() {

    log.info"""
    Usage:
      nextflow run ./crisprSNF --input design.csv 
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

ch_multiqc_config        = ''
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

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

params.publish_dir_mode = 'copy'

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
    ch_samtools_stats    = DEDUP_UMI_UMITOOLS.out.stats
    ch_samtools_flagstat = DEDUP_UMI_UMITOOLS.out.flagstat
    ch_samtools_idxstats = DEDUP_UMI_UMITOOLS.out.idxstats

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
    ch_grna_counts = COUNT_GRNA_READS.out.counts
    ch_versions = ch_versions.mix(COUNT_GRNA_READS.out.versions)

    // 07 From here we have a table 

     //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}


workflow {
    CRISPRSNF ()
}
