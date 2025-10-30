#!/usr/bin/env nextflow

process PANGENIE {

    container 'svcalling-pipeline:latest'
    
    tag "pangenie.${sample_id}"
    //publishDir "${params.output_dir}/${sample_id}/variants", mode: 'copy'
    stageInMode 'symlink'

    cpus 24
    memory '50GB'

    input:
    path pangenome_dir
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val (sample_id), path("${sample_id}_result_SV_genotyping.vcf"), emit: pangenie_gt_calls

    script:
    /*
    -j is used to define the number of threads for kmer-counting (default = 1)
    -t is used to define the number of threads for the core algorithm (default = 1, max = nChromosomes in vcf)
    */
    def third = task.cpus / 3
    println "Using $third cpus for kmer counting..."
    println "Using $task.cpus for core algorithm..."
    """
    PanGenie -f $params.pangenome_data_prefix \
    -i <(zcat ${fastq1} ${fastq2}) \
    -o ${sample_id}_result_SV \
    -s $sample_id \
    -j $third \
    -t $task.cpus 
    """

}