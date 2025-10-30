#!/usr/bin/env nextflow


// Process: Generate metrics
process METRICS {
    tag "metrics.${sample_id}"
    publishDir "${dest_prefix}/${sample_id}", mode: 'copy'

    cpus 2
    memory '8GB'

    input:
    tuple val(sample_id), path(bam), path(bai)
    tuple path(fasta), path(indexes)
    path(target_bed)
    val dest_prefix

    output:
    tuple val(sample_id), path("${sample_id}_metrics.tar.gz"), path("${sample_id}_QCMetrics.json"), emit: metrics_qc

    script:
    /*
    region.interval_list for WGS analysis is canonical protien coding gene regions from gencode v49
    
    ${params.sentieon_dir}/bin/sentieon plot GCBias -o gc-report.pdf gc_metrics.txt
    ${params.sentieon_dir}/bin/sentieon plot QualDistribution -o qd-report.pdf qd_metrics.txt
    ${params.sentieon_dir}/bin/sentieon plot MeanQualityByCycle -o mq-report.pdf mq_metrics.txt
    ${params.sentieon_dir}/bin/sentieon plot InsertSizeMetricAlgo -o is-report.pdf is_metrics.txt
    */
    """
    echo -e '@HD\tVN:1.0\tSO:coordinate' > region.interval_list
    ${params.samtools_dir} view -H ${bam} | grep '^@SQ' >> region.interval_list
    grep -v ^# ${target_bed} | awk -F'\t' 'BEGIN{OFS=FS}{print \$1,\$2+1,\$3,"+",\$4;}' >> region.interval_list

    ${params.sentieon_dir}/bin/sentieon driver -r ${fasta} -t ${task.cpus} -i ${bam} \
        --algo MeanQualityByCycle mq_metrics.txt \
        --algo QualDistribution qd_metrics.txt \
        --algo GCBias --summary gc_summary.txt gc_metrics.txt \
        --algo AlignmentStat --adapter_seq '' aln_metrics.txt \
        --algo InsertSizeMetricAlgo is_metrics.txt \
        --algo HsMetricAlgo --targets_list region.interval_list --baits_list region.interval_list hs_metrics.txt
    
    sentieon-metrics2json -a aln_metrics.txt \
        -g gc_metrics.txt \
        -u gc_summary.txt \
        -h hs_metrics.txt \
        -i is_metrics.txt \
        -m mq_metrics.txt \
        -q qd_metrics.txt > ${sample_id}_QCMetrics.json

    mkdir -p metrics
    mv gc_metrics.txt gc_summary.txt aln_metrics.txt qd_metrics.txt mq_metrics.txt is_metrics.txt hs_metrics.txt metrics/
    tar -czf ${sample_id}_metrics.tar.gz metrics
    """
}

