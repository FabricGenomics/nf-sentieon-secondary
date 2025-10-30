#!/usr/bin/env nextflow

process ADD_TRF {

    container 'truvari:latest'
    
    tag "TRFannotation.${sample_id}"
    stageInMode 'symlink'
    //publishDir "${params.output_dir}/${sample_id}/variants", mode: 'copy'

    cpus 4
    memory '24GB'

    input:
    tuple path(fasta), path(indexes)
    tuple val(sample_id), path(input_vcf), path(input_vcf_idx)

    output:
    tuple val(sample_id), path("${sample_id}_TRannotated_snv_indel_sv.vcf.gz"), path("${sample_id}_TRannotated_snv_indel_sv.vcf.gz.tbi"), emit: final_sv_calls

    script:
    """
    truvari anno trf -i $input_vcf \
     -o ${sample_id}_TRannotated.vcf \
     -e /opt/truvari-source/trf409.linux64 \
     -f $fasta \
     -r /opt/truvari-source/anno.trf.bed.gz \
     -t $task.cpus 

     bcftools sort ${sample_id}_TRannotated.vcf -o ${sample_id}_TRannotated_snv_indel_sv.vcf

     bgzip ${sample_id}_TRannotated_snv_indel_sv.vcf
     tabix ${sample_id}_TRannotated_snv_indel_sv.vcf.gz

     # Remove intermediate files
     rm ${sample_id}_TRannotated.vcf 
    """

}