#!/usr/bin/env nextflow

process MANTA_HS {

    container 'manta-hs:latest'
    
    tag "manta_SVs.${sample_id}"
    stageInMode 'symlink'
    cpus 8
    memory '24GB'

    input:
    tuple path(fasta), path(indexes)
    tuple val(sample_id), path(bam), path(bai)

    output:
    path("${sample_id}_manta_out/results/variants/diploidSV.vcf"), emit: manta_vcf

    script:
    """
    /usr/local/bin/manta/bin/configManta.py \
    --bam $bam \
    --referenceFasta $fasta \
    --runDir ${sample_id}_manta_out

    ${sample_id}_manta_out/runWorkflow.py -j $task.cpus
    gunzip ${sample_id}_manta_out/results/variants/diploidSV.vcf.gz
    """

}