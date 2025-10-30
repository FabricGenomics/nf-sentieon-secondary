#!/usr/bin/env nextflow

// Process: Mapping reads with BWA-MEM, sorting, and indexing
process BWA_MEM {
    tag "mem.${sample_id}" // labels log line as mem.sampleA instead of (1)

    cpus 45
    memory '64GB'

    //publishDir "${params.output_dir}/${sample_id}", mode: 'copy', pattern: "*.{bam,bai}" // copies output to publishDir
    stageInMode 'symlink' // copies input files into work directory 

    input:
      tuple val(sample_id), path(fastq_1), path(fastq_2)
      path dnascope_bundle
      tuple path(fasta), path(indexes)


    output:
      tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: sorted_bam
    // tuple val(sample_id), path("${sample_id}.sam"), emit: aligned_sam


    script:
    """
    ${params.sentieon_dir}/bin/sentieon bwa mem \
        -R "@RG\\tID:rg_${sample_id}\\tSM:${sample_id}\\tPL:${params.platform}" \
        -t ${task.cpus} -K 10000000 -x ${dnascope_bundle}/${params.bwa_model} \
        ${fasta} ${fastq_1} ${fastq_2} |\
        ${params.sentieon_dir}/bin/sentieon util sort -i - -r ${fasta} \
        -o ${sample_id}_sorted.bam -t ${task.cpus} --temp_dir /home/ec2-user/tmp --sam2bam

      # Index the BAM file
    ${params.sentieon_dir}/bin/sentieon util index ${sample_id}_sorted.bam

    """

}