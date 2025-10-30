#!/usr/bin/env nextflow

process CONCAT_FINAL {

    container 'svcalling-pipeline:latest'
    
    tag "concatenate_vcfs.${sample_id}"
    publishDir "${dest_prefix}/${sample_id}", mode: 'copy'
    stageInMode 'symlink'

    cpus 1
    memory '16GB'

    input:
    tuple val (sample_id), path(concat1_vcf), path(concat1_vcf_tbi)
    tuple val (sample_id_copy), path(cnv_vcf), path(cnv_vcf_tbi)
    val dest_prefix


    output:
    tuple val(sample_id), path("${sample_id}_allVars_concat.vcf.gz"), path("${sample_id}_allVars_concat.vcf.gz.tbi"), emit: final_concat_vcf
    
    script:
    """
    bcftools concat -d exact \
    -a -Oz -o ${sample_id}_allVars_concat.vcf.gz \
    $concat1_vcf $cnv_vcf 

    tabix ${sample_id}_allVars_concat.vcf.gz
    """
}