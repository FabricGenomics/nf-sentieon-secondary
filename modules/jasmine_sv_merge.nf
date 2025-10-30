#!/usr/bin/env nextflow

process JASMINE_SV_MERGE {

    container 'jasmine:latest'
    containerOptions '-u 1000:1000'
    
    tag "merge_SVs.${sample_id}"
    stageInMode 'symlink'
    cpus 8
    memory '24GB'

    input:
    tuple path(fasta), path(indexes)
    tuple val(sample_id), path(pangenie_vcf)
    path(manta_vcf)
    path(delly_vcf)

    output:
    tuple val(sample_id), path("${sample_id}_jasmine_merge.vcf"), emit: jasmine_out_vcf

    script:
    """
    ls $delly_vcf >> files.txt
    ls $manta_vcf >> files.txt
    ls $pangenie_vcf >> files.txt

    # Call Jasmine
    jasmine file_list=files.txt \
    out_file=${sample_id}_jasmine_merge.vcf \
    genome_file=$fasta \
    --allow_intrasample \
    --output_genotypes \
    --ignore_strand \
    --dup_to_ins \
    --centroid_merging \
    threads=8
    """

}