#!/usr/bin/env nextflow

process PROCESS_PANGENIE {

    container 'svcalling-pipeline:latest'
    
    tag "pangenie_processing.${sample_id}"
    stageInMode 'symlink'

    cpus 8
    memory '50GB'

    input:
    tuple path(ref_vcf), path(ref_vcf_idx)
    tuple val (sample_id), path(input_vcf)

    output:
    tuple val(sample_id), path("${sample_id}_pangenie_jasmine_input.vcf"), emit: pangenie_input_jasmine
    //tuple path("${sample_id}_pangenie_future_merge.vcf.gz"), path("${sample_id}_pangenie_future_merge.vcf.gz.tbi"), emit: pangenie_final_merge
    script:
    /*
    # Sort, filter and compress all Pangenie calls for truvari
    bcftools view -i '(GT != "RR" && GT != "mis")' \
    -Oz -o ${sample_id}_pangenie_future_merge.vcf.gz \
    ${sample_id}_result_SV_genotyping_biallelic_SVinfo_sorted.vcf.gz

    tabix ${sample_id}_pangenie_future_merge.vcf.gz
    */
    """
    # Create biallelic VCF output from pangenie - script inside pangenie build in container
    cat $input_vcf |\
    python3 /repos/pangenie/pipelines/run-from-callset/scripts/convert-to-biallelic.py $ref_vcf > ${sample_id}_result_SV_genotyping_biallelic.vcf
    
    # Add SV INFO fields to VCF
    sv_info_field.py \
    --input ${sample_id}_result_SV_genotyping_biallelic.vcf \
    --output ${sample_id}_result_SV_genotyping_biallelic_SVinfo.vcf

    bgzip ${sample_id}_result_SV_genotyping_biallelic_SVinfo.vcf
    tabix ${sample_id}_result_SV_genotyping_biallelic_SVinfo.vcf.gz

    # Sort 
    bcftools sort ${sample_id}_result_SV_genotyping_biallelic_SVinfo.vcf.gz \
    -o ${sample_id}_result_SV_genotyping_biallelic_SVinfo_sorted.vcf.gz

    tabix ${sample_id}_result_SV_genotyping_biallelic_SVinfo_sorted.vcf.gz

    # Filter for jasmine input
    bcftools view -i '(SVLEN < -49 || SVLEN > 49) && (GT != "RR" && GT != "mis")' \
    -Ov -o ${sample_id}_pangenie_jasmine_input.vcf \
    ${sample_id}_result_SV_genotyping_biallelic_SVinfo_sorted.vcf.gz

    # rm intermediate files
    rm ${sample_id}_result_SV_genotyping_biallelic.vcf 
    rm ${sample_id}_result_SV_genotyping_biallelic_SVinfo.vcf.gz ${sample_id}_result_SV_genotyping_biallelic_SVinfo.vcf.gz.tbi
    rm ${sample_id}_result_SV_genotyping_biallelic_SVinfo_sorted.vcf.gz ${sample_id}_result_SV_genotyping_biallelic_SVinfo_sorted.vcf.gz.tbi
    """

}