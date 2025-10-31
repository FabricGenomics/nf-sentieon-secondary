#!/usr/bin/env nextflow
// Process: CNV calling

process CNV_CALLING_WES {
    tag "cnv_calling.${sample_id}"
    //publishDir "${params.output_dir}/${sample_id}/variants", mode: 'copy'

    cpus 2
    memory '8GB'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path sentieon_bundle
    tuple path(fasta), path(indexes)
    path target_bed

    output:
    tuple val(sample_id), path("${sample_id}_CNV.vcf.gz"), path("${sample_id}_CNV.vcf.gz.tbi"), emit: cnv_vcf

    script:
    """
    ${params.sentieon_dir}/bin/sentieon driver -t ${task.cpus} -r ${fasta} -i ${bam} \
        --interval ${target_bed} --algo CNVscope \
        --model ${sentieon_bundle}/${params.cnv_model} tmp_cnv.vcf.gz

    ${params.sentieon_dir}/bin/sentieon driver -t ${task.cpus} -r ${fasta} \
        --algo CNVModelApply --model ${sentieon_bundle}/${params.cnv_model} \
        -v tmp_cnv.vcf.gz ${sample_id}_CNV.vcf.gz
    """
}
