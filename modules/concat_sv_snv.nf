#!/usr/bin/env nextflow

process CONCAT1 {

    container 'svcalling-pipeline:latest'
    
    tag "concat_SNV_SV.${sample_id}"
    stageInMode 'symlink'

    cpus 1
    memory '16GB'

    input:
    tuple val (sample_id), path(sv_vcf), path(sv_vcf_tbi)
    tuple val (sample_id_copy), path(smallvar_vcf), path(smallvar_vcf_tbi)

    output:
    tuple val(sample_id), path("${sample_id}_snv_indel_sv.vcf.gz"), path("${sample_id}_snv_indel_sv.vcf.gz.tbi"), emit: concat_vcf
    
    script:
    """
    tabix $smallvar_vcf

    bcftools concat -d exact \
    -a -Oz -o ${sample_id}_snv_indel_sv.vcf.gz \
    $sv_vcf $smallvar_vcf 

    tabix ${sample_id}_snv_indel_sv.vcf.gz
    """
}