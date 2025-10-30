#!/usr/bin/env nextflow

process CNV_PROCESSING {

    container 'svcalling-pipeline:latest'
    
    tag "cnv_processing.${sample_id}"
    //publishDir "${params.output_dir}/${sample_id}/variants", mode: 'copy'
    stageInMode 'symlink'
    cpus 1
    memory '8GB'

    input:
    tuple val(sample_id), path(input_vcf)

    output:
    tuple val(sample_id), path("${sample_id}_CNV.vcf.gz"), path("${sample_id}_CNV.vcf.gz.tbi"), emit: cnv_vcf

    script:
    """
    echo "${sample_id}_cnv" >> samps.txt
    echo "${sample_id}_sv" >> samps.txt
    echo "${sample_id}" >> final_samp.txt

    cnv_vcf_update.py -i ${input_vcf} -o ${sample_id}_CNV_processed.vcf

    bcftools sort ${sample_id}_CNV_processed.vcf -o ${sample_id}_CNV_processed_sorted.vcf
    bcftools reheader -s samps.txt -o ${sample_id}_CNV_multisamp.vcf ${sample_id}_CNV_processed_sorted.vcf

    bgzip ${sample_id}_CNV_multisamp.vcf  
    tabix ${sample_id}_CNV_multisamp.vcf.gz

    bcftools view -s ${sample_id}_cnv ${sample_id}_CNV_multisamp.vcf.gz -Oz |\
    bcftools reheader -s final_samp.txt -o ${sample_id}_CNV.vcf.gz 

    tabix ${sample_id}_CNV.vcf.gz

    rm ${sample_id}_CNV_processed.vcf ${sample_id}_CNV_processed_sorted.vcf 
    rm samps.txt final_samp.txt
    rm ${sample_id}_CNV_multisamp.vcf.gz ${sample_id}_CNV_multisamp.vcf.gz.tbi
    """

}