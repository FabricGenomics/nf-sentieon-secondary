#!/usr/bin/env nextflow

process JASMINE_CNV {

    container 'jasmine:latest'
    containerOptions '-u 1000:1000'
    
    tag "merge_CNVs.${sample_id}"
    stageInMode 'symlink'
    cpus 8
    memory '16GB'

    input:
    tuple path(fasta), path(indexes)
    path(delly_vcf)
    tuple val(sample_id), path(delly_cnv)

    output:
    tuple val(sample_id), path("${sample_id}_jasmine_merge_cnv.vcf"), emit: jasmine_cnv_out_vcf

    script:
    """
    ls ${delly_cnv} >> files.txt
    ls ${delly_vcf} >> files.txt
    
    # Call Jasmine
    jasmine file_list=files.txt \
    out_file=${sample_id}_jasmine_merge_cnv.vcf \
    genome_file=$fasta \
    max_dist=20 \
    min_support=2 \
    threads=8 \
    --use_end \
    --output_genotypes \
    --ignore_type \
    --ignore_strand
    """

}