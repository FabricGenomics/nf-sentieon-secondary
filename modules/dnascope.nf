#!/usr/bin/env nextflow

// Process: DNAscope variant calling
process DNASCOPE {
    tag "variant_calling.${sample_id}"
    //publishDir "${params.output_dir}/${sample_id}/variants", mode: 'copy'

    cpus 48 
    memory '16GB'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path dnascope_bundle
    tuple path(fasta), path(indexes)

    output:
    tuple val(sample_id), path("${sample_id}_output-ds.vcf.gz"), emit: dnascope_vcf

    script:
    def indel_model = params.pcr_free == "true" ? "--pcr_indel_model none" : ""
    """
    ${params.sentieon_dir}/bin/sentieon driver -r ${fasta} -t ${task.cpus} -i ${bam} \
        --algo DNAscope ${indel_model} --model ${dnascope_bundle}/${params.dnascope_model} \
        output-ds_tmp.vcf.gz

    ${params.sentieon_dir}/bin/sentieon driver -r ${fasta} -t ${task.cpus} \
        --algo DNAModelApply --model ${dnascope_bundle}/${params.dnascope_model} \
        -v output-ds_tmp.vcf.gz ${sample_id}_output-ds.vcf.gz
    """
}