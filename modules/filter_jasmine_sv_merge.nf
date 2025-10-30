#!/usr/bin/env nextflow

process FILTER_JASMINE_SVMERGE {

    container 'svcalling-pipeline:latest'
    
    tag "filter_mergedSVs.${sample_id}"
    //publishDir "${params.output_dir}/${sample_id}/variants", mode: 'copy'
    stageInMode 'symlink'
    cpus 4
    memory '16GB'

    input:
    tuple val(sample_id), path(jasmine_vcf)
    tuple val(sample_id_cp), path(pangenie_vcf)
    path header_lines

    output:
    tuple val(sample_id), path("${sample_id}_combined_svs.vcf.gz"), path("${sample_id}_combined_svs.vcf.gz.tbi"), emit: combined_svs

    script:
    """
    # Create new sample header line with same order as read into jasmine 
    echo "delly" >> sample_names.txt
    echo "manta" >> sample_names.txt
    echo "pangenie" >> sample_names.txt
    
    # Create new sample header line for eventual subset file 
    echo "${sample_id}" > "${sample_id}".txt

    # Correct header to contain all INFO and FORMAT fields
    bcftools annotate -h $header_lines $jasmine_vcf |\
    bcftools reheader -s sample_names.txt |\
    bcftools sort -Oz -o jasmine_reheadered_sorted.vcf.gz

    tabix jasmine_reheadered_sorted.vcf.gz
    
    # Identify SVs absent from PanGenome
    filter_vcf_samps.py -i jasmine_reheadered_sorted.vcf.gz\
    -o jasmine_noPangenie_calls.vcf --sample-a pangenie --sample-b manta --sample-c delly

    # Subset to Delly GTs only and rename sample 
    bcftools view -s delly -Oz -o delly_only_gts.vcf.gz jasmine_noPangenie_calls.vcf
    tabix delly_only_gts.vcf.gz

    bcftools reheader -s "${sample_id}".txt -o delly_only_gts_${sample_id}.vcf.gz delly_only_gts.vcf.gz
    tabix delly_only_gts_${sample_id}.vcf.gz

    # Combine with all PanGenie calls
    bcftools view -Oz -o pangenie_svs.vcf.gz $pangenie_vcf
    tabix pangenie_svs.vcf.gz

    bcftools concat -a delly_only_gts_${sample_id}.vcf.gz pangenie_svs.vcf.gz |\
    bcftools sort -Oz -o ${sample_id}_combined_svs.vcf.gz
    tabix ${sample_id}_combined_svs.vcf.gz

    # Remove tmp files
    rm jasmine_reheadered_sorted.vcf.gz \
    jasmine_reheadered_sorted.vcf.gz.tbi \
    jasmine_noPangenie_calls.vcf \
    delly_only_gts.vcf.gz \
    delly_only_gts.vcf.gz.tbi \
    delly_only_gts_${sample_id}.vcf.gz \
    delly_only_gts_${sample_id}.vcf.gz.tbi \
    pangenie_svs.vcf.gz \
    pangenie_svs.vcf.gz.tbi
    """
}