#!/usr/bin/env nextflow

process DELLY_SV {

    container 'delly:latest'
    
    tag "delly_SVs.${sample_id}"
    stageInMode 'symlink'
    cpus 2
    memory '48GB'

    input:
    tuple path(fasta), path(indexes)
    tuple val(sample_id), path(bam), path(bai)

    output:
    path("${sample_id}_delly.vcf"), emit: delly_vcf

    script:
    """
    delly call -g $fasta $bam -x /maps/hg38_exclude.tsv -o ${sample_id}_delly.bcf
    bcftools view -Ov -o ${sample_id}_delly.vcf ${sample_id}_delly.bcf
    rm ${sample_id}_delly.bcf
    """

}