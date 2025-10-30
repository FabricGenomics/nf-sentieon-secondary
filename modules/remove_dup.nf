#!/usr/bin/env nextflow

// Process: Remove duplicates
process DEDUP {
    tag "dedup.${sample_id}"
    publishDir "${dest_prefix}/${sample_id}", mode: 'copy', pattern: "*.{bam,bai}"

    cpus 16
    memory '32GB'

    input:
    tuple val(sample_id), path(bam), path(bai)
    val dest_prefix

    output:
    tuple val(sample_id), path("${sample_id}_deduped.bam"), path("${sample_id}_deduped.bam.bai"), emit: deduped_bam
    path "${sample_id}_dedup_metrics.txt", emit: dedup_metrics

    script:
    """
    ${params.sentieon_dir}/bin/sentieon driver -t ${task.cpus} -i ${bam} \
        --algo LocusCollector --fun score_info score.txt

    ${params.sentieon_dir}/bin/sentieon driver -t ${task.cpus} -i ${bam} \
        --algo Dedup --score_info score.txt \
        --metrics ${sample_id}_dedup_metrics.txt ${sample_id}_deduped.bam

    # Index the deduped BAM file
    ${params.sentieon_dir}/bin/sentieon util index ${sample_id}_deduped.bam
    """
}