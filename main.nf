#!/usr/bin/env nextflow

// Enable DSL2 - this is the domain specific language in nextflow, hello-nextflow used this
nextflow.enable.dsl=2

// Parameters from command line - default is null
params.proband_fastq1 = null 
params.proband_fastq2 = null
params.mother_fastq1 = null
params.mother_fastq2 = null
params.father_fastq1 = null
params.father_fastq2 = null
params.s3_bucket_out = params.s3_bucket_out
    ? params.s3_bucket_out.replaceAll(/\/+$/, '')
    : null
params.wes = false
params.platform = null
// Config file parameters
params.dnascope_model = null
params.bwa_model = null
params.cnv_model = null
params.sentieon_complete_mgi_wgs_bundle = null
params.sentieon_complete_mgi_wes_bundle = null
params.sentieon_element_wgs_bundle = null
params.sentieon_element_wes_bundle = null
params.sentieon_illumina_wgs_bundle = null
params.sentieon_illumina_wes_bundle = null
params.pcr_free = null
params.reference_fasta = null
params.sentieon_dir = null
params.samtools_dir = null
params.pangenome_data_dir = null
params.pangenome_data_prefix = null
params.biallelic_vcf = null
params.target_bed = null

// import processes from './modules/'
include { BWA_MEM } from './modules/bwamem_process.nf'
include { METRICS_WES } from './modules/metrics_wes.nf'
include { DEDUP } from './modules/remove_dup.nf'
include { DNASCOPE } from './modules/dnascope.nf'
include {DNASCOPE_WES } from './modules/dnascope_wes.nf'
include { PANGENIE } from './modules/pangenie.nf'
include { PROCESS_PANGENIE } from './modules/pangenie_processing.nf'
include { ADD_TRF } from './modules/trf_anno.nf'
include { DELLY_SV } from './modules/delly_sv.nf'
include { DELLY_CNV } from './modules/delly_cnv.nf'
include { JASMINE_CNV } from './modules/jasmine_cnv.nf'
include { GANGSTR } from './modules/gangstr.nf'
include { CNV_PROCESSING } from './modules/cnv_processing.nf'
include { CNV_CALLING_WES } from './modules/cnv_calling_wes.nf'
include { CONCAT1 } from './modules/concat_sv_snv.nf'
include { CONCAT_FINAL } from './modules/concat_final.nf'

// Workflow
workflow {
    // Validate required parameters
    if (!params.proband_fastq1 || !params.proband_fastq2) {
        error "Error: Proband FASTQ files (--proband_fastq1 and --proband_fastq2) are required."
    }
    if (!params.s3_bucket_out) {
        error "Error: Secondary S3 Bucket (--s3_bucket_out) is required."
    }
    if (!params.reference_fasta) { error "Error: Please specify reference fasta with --reference_fasta" }
    if (!params.bwa_model) { error "Error: Please specify bwa model with --bwa_model" }
    if (!params.sentieon_dir) { error "Error: Please specify sentieon dir with --sentieon_dir" }
    if (!params.samtools_dir) { error "Error: Please specify samtools dir with --samtools_dir" }
    if (!params.pcr_free) { error "Error: Please specify pcr free with --pcr_free true / false" }
    if (!params.platform) { error "Error: Sequencing platform is required. Please specify platform with --platform" }

    log.info "Workflow session ID: ${workflow.sessionId}"

    // Generate expected input for couchDB
    def runId = workflow.sessionId

    def bucketPrefix = (params.s3_bucket_out ?: '').replaceAll(/\/+$/, '')
    log.info "bucketPrefix = ${bucketPrefix}"
    def destPrefix = "${bucketPrefix}/${runId}"
    dest_prefix_ch = Channel.value(destPrefix)

    // Define sample IDs
    def proband_id = "proband"
    def mother_id = params.mother_fastq1 ? "mother" : null
    def father_id = params.father_fastq1 ? "father" : null

    // Create channels for input FASTQ files
    input_fastqs = Channel.from([[
        id: proband_id,
        fastq1: file(params.proband_fastq1),
        fastq2: file(params.proband_fastq2)
    ]])
    if (mother_id) {
        input_fastqs = input_fastqs.mix(Channel.from([[
            id: mother_id,
            fastq1: file(params.mother_fastq1),
            fastq2: file(params.mother_fastq2)
        ]]))
    }
    if (father_id) {
        input_fastqs = input_fastqs.mix(Channel.from([[
            id: father_id,
            fastq1: file(params.father_fastq1),
            fastq2: file(params.father_fastq2)
        ]]))
    }

    // Bundle inputs
    if (params.platform == "CG" || params.platform == "MGI") {
        dnascope_bundle_ch = Channel.value(params.sentieon_complete_mgi_wgs_bundle)
        wes_dnascope_bundle_ch = Channel.value(params.sentieon_complete_mgi_wes_bundle)
    } else if (params.platform == "Element") {
        dnascope_bundle_ch = Channel.value(params.sentieon_element_wgs_bundle)
        wes_dnascope_bundle_ch = Channel.value(params.sentieon_element_wes_bundle)
    } else if (params.platform == "Illumina") {
        dnascope_bundle_ch = Channel.value(params.sentieon_illumina_wgs_bundle)
        wes_dnascope_bundle_ch = Channel.value(params.sentieon_illumina_wes_bundle)
    } else {
        error "Unknown Sequencing Platform Type: ${params.platform}"
    }

    // FASTAs and indexes
    def index_exts = ['fai','amb','ann','bwt','pac','sa']
    def fasta_indexes = index_exts.collect { ext -> file("${params.reference_fasta}.${ext}") }
    reference_fasta_ch = Channel.value( tuple(file(params.reference_fasta), fasta_indexes) )

    // STR BED file 
    gangstr_bed_ch = Channel.fromPath('resources/hg38_ver13.bed', checkIfExists: true)

    // Run pipeline based on WES or WGS input data
    if (params.wes) {
        println "Running pipeline for WES input."

        if (!params.target_bed) { 
            error "Error: Bed file containing target regions (--target_bed) is required for WES mode."
        }
        if (!params.cnv_model) { 
            error "Error: Sentieon CNV model (--cnv_model) is required for WES mode."
        }

        //target bed file channel
        target_bed_ch = Channel.fromPath(params.target_bed, checkIfExists: true)

        // Map FASTQ files to BAM
        bwa_mem_out = BWA_MEM(input_fastqs.map { [it.id, it.fastq1, it.fastq2] } , wes_dnascope_bundle_ch, reference_fasta_ch )
        bwa_mem_out.sorted_bam.view()
        
        // WES Metrics
        metrics_out = METRICS_WES(bwa_mem_out.sorted_bam, reference_fasta_ch, target_bed_ch, dest_prefix_ch)

        // Remove duplicates
        dedup_out = DEDUP(bwa_mem_out.sorted_bam, dest_prefix_ch)

        // DNAscope WES
        dnascope_wes_out = DNASCOPE_WES(dedup_out.deduped_bam, wes_dnascope_bundle_ch, reference_fasta_ch, target_bed_ch)

        // CNV Calling WES
        cnv_wes_out = CNV_CALLING_WES(dedup_out.deduped_bam, dnascope_bundle_ch, reference_fasta_ch, target_bed_ch)

        // Combine VCF
        concat_out_wes = CONCAT_FINAL(dnascope_wes_out.dnascope_vcf, cnv_wes_out.cnv_vcf, dest_prefix_ch)

    } else {
        println "Running pipeline for WGS input."
        
        if (!params.pangenome_data_dir) { 
            error "Error: Directory containing PanGenome input files (--pangenome_data_dir) is required for WGS mode." 
        }
        if (!params.pangenome_data_prefix) { 
            error "Error: The PanGenome data prefix (--pangenome_data_prefix) is required for WGS mode."
        }
        if (!params.biallelic_vcf) { 
            error "Error: Biallelic PanGenome VCF (--biallelic_vcf) is required for WGS mode."
        }

        // Pangenie Input
        pangenie_ref_ch = Channel.fromPath("${params.pangenome_data_dir}/*").collect()

        // Pangenie Output - add SV Info 
        pangenie_py_ch = Channel.of(
            tuple(
                file(params.biallelic_vcf),
                file("${params.biallelic_vcf}.csi")
                )
            )
        // Map FASTQ files to BAM
        bwa_mem_out = BWA_MEM(input_fastqs.map { [it.id, it.fastq1, it.fastq2] } , dnascope_bundle_ch, reference_fasta_ch )
        bwa_mem_out.sorted_bam.view()
        
        // Remove duplicates
        dedup_out = DEDUP(bwa_mem_out.sorted_bam, dest_prefix_ch)

        // DNAscope variant calling
        dnascope_out = DNASCOPE(dedup_out.deduped_bam, dnascope_bundle_ch, reference_fasta_ch)

        // Structural variant genotyping with pangenie
        pangenie_out = PANGENIE(pangenie_ref_ch, input_fastqs.map { [it.id, it.fastq1, it.fastq2] })
        pangenie_processing_out = PROCESS_PANGENIE(pangenie_py_ch, pangenie_out.pangenie_gt_calls)

        // Call STRs with GangSTR
        str_calls = GANGSTR(dedup_out.deduped_bam, reference_fasta_ch, gangstr_bed_ch)
        
        // Combine SNV/INDEL VCF and SV VCF for STR/VNTR annotation
        concat_sv_smallvar = CONCAT1(pangenie_processing_out.filtered_SVs, dnascope_out.dnascope_vcf)

        // Run Truvari TRF annotation of SNVs, INDELs, and SVs 
        trf_annotation_out = ADD_TRF(reference_fasta_ch, concat_sv_smallvar.concat_vcf)

        // CNV calling
        delly_out = DELLY_SV(reference_fasta_ch, dedup_out.deduped_bam)
        cnv_out = DELLY_CNV(reference_fasta_ch, dedup_out.deduped_bam, delly_out.delly_vcf)
        merged_out = JASMINE_CNV(reference_fasta_ch, delly_out.delly_vcf, cnv_out.delly_cnv_vcf)
        cnv_calls = CNV_PROCESSING(merged_out.jasmine_cnv_out_vcf)

        //Final concatenation step
        final_vcf = CONCAT_FINAL(trf_annotation_out.final_sv_calls, cnv_calls.cnv_vcf, dest_prefix_ch)
    }
}
