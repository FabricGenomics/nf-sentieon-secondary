#!/usr/bin/env nextflow

process GANGSTR {

    container 'gymreklab/str-toolkit:latest'
    
    tag "STRcalling.${sample_id}"
    //publishDir "${params.output_dir}/${sample_id}/variants", mode: 'copy'
    stageInMode 'symlink'

    cpus 16
    memory '32GB'

    input:
    tuple val(sample_id), path(bam), path(bai)
    tuple path(fasta), path(indexes)
    path str_bed

    output:
    tuple val (sample_id), path("${sample_id}.gangSTR.dumpSTR.vcf"), emit: gangstr_vcf

    script:
    """
    GangSTR --bam ${bam} --ref ${fasta} \
    --regions ${str_bed} --out ${sample_id}.gangSTR

    dumpSTR --vcf ${sample_id}.gangSTR.vcf \
    --out ${sample_id}.gangSTR.dumpSTR \
    --filter-spanbound-only --filter-badCI
    """

}