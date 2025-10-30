#!/usr/bin/env nextflow

process DELLY_CNV {

    container 'delly:latest'
    
    tag "delly_CNVs.${sample_id}"
    stageInMode 'symlink'
    cpus 2
    memory '24GB'

    input:
    tuple path(fasta), path(indexes)
    tuple val(sample_id), path(bam), path(bai)
    path(sv_vcf)

    output:
    tuple val(sample_id), path("${sample_id}_delly_cnv.vcf"), emit: delly_cnv_vcf

    script:
    /*
    Should we include bed file of gene regions to target for CNV calling?
    */
    """
    delly cnv -g ${fasta} \
    -m /maps/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz \
    -i 1000 \
    -j 1000 \
    -l ${sv_vcf} \
    -o ${sample_id}_delly_cnv.bcf ${bam}

    bcftools annotate -I "%ID\\_copy_num_[%CN]" -Ov -o ${sample_id}_delly_cnv.vcf ${sample_id}_delly_cnv.bcf
    #rm ${sample_id}_delly_cnv.bcf
    """

}
