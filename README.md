# feature/wes_update

## Concept 
The updated Secondary Analysis pipeline can handle WGS or WES data. For both WGS and WES, the Sentieon driver algorithms are used for calling SNVs and INDELs. For WGS, we use a PanGenome approach that utilizes PanGenie to genotype SVs that are present in a PanGenome containing 88 haplotypes. PanGenie was found to increase recall to ~91% from 15%-30% depending on the caller. For CNVs, delly SV + delly CNV callers are used. For STRs, the PanGenie+SentieonDNAscope calls are annotated with `truvari anno trf` which uses Tandem Repeat Finder to annotate tandem repeats in the variant calls. For WES, we do not call SVs. CNVs are called with Sentieon's CNVscope algorithm, but it is recommended to use WGS for CNVs if possible. 

Improvement of the SV/CNV calling approach is linked to the JIRA epic [FE3-2205](https://spiralgenetics.atlassian.net/browse/FE3-2205). 

## Data Locations
**Benchmarking**
* GiaB HPRC Draft SV Benchmark: "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/" 
* GiaB HPRC medical gene SV BenchmarK: "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/"  

**HG002 Sequencing Data** 
* HG002 sequencing data for testing: "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/" 
* HG002 Illumina fastq files: "gs://brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x/"  

**SV Caller**
* Pangenie: "https://github.com/eblerjana/pangenie/tree/master"  

**Reference Files** 
* Human Genome (GRCh38): "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"  
* PanGenome reference VCF (88 haplotypes): "https://zenodo.org/record/6797328/files/cactus_filtered_ids.vcf.gz"
* PanGenome reference biallelic-callset: "https://zenodo.org/record/6797328/files/cactus_filtered_ids_biallelic.vcf.gz"  
* Gencode v49 (GRCh38.p14) GFF3 file used for creating the protein-coding gene regions target bed file: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gff3.gz"  

**Other files**
* Tandem Repeat Finder executable: "https://github.com/Benson-Genomics-Lab/TRF/releases/tag/v4.09.1"  
* Adotto TR Catalog: "https://github.com/ACEnglish/adotto/blob/main/regions/DataDescription.md"  
* Sentieon DNAscope bundles: "https://github.com/Sentieon/sentieon-models?tab=readme-ov-file" 

## Description of Files. 
* modules/
    - bwamem_process.nf aligns FASTQ files with BWA 
    - cnv_processing.nf filters the jasmine merged CNV calls and produces the final CNV VCF 
    - concat_final.nf produces the final concatenated VCF file to be read into FE3 
    - concat_sv_snv.nf concatenates the SNV/INDEL and PanGenie SV VCF for tandem repeat annotations 
    - delly_cnv.nf calls CNVs with Delly and uses Delly SV calls to refine breakpoints 
    - delly_sv.nf calls SVs with Delly to be used in CNV calling and merged with Manta and Pangenie SV calls 
    - dnascope.nf calls sentieon for SNV/INDEL calling 
    - jasmine_cnv.nf merges Delly CNV and Delly SV calls to obtain DEL/DUP information 
    - metrics.nf outputs the expected metrics tarball and JSON files expected by the API 
    - remove_dup.nf marksand removes duplicate reads 
    - pangenie.nf calls pangenie to genotype SVs against the pangenome 
    - pangenie_processing.nf converts the SV VCF to a biallelic format and adds SV LENGTH, SV TYPE, and END fields to INFO. 
    - trf_anno.nf runs is the nextflow process that runs `truvari anno trf` to find and annotate tandem repeats 
    - cnv_calling_wes.nf is the nextflow process that calls CNVs for WES data  
    - dnascope_wes.nf is the nextflow process that calls SNVs/INDELs for WES data  
    - metrics_wes.nf calculates metrics for WES data including HS metrics  
* bin/ 
    - cnv_vcf_update.py is a python script to add copy number info to the INFO field and update the SVTYPE to DEL or DUP for CNV calls 
    - sv_info_field.py is a custom python script to add common SV annotations to the VCF output by pangenie and is called by nextflow 
    - extract_trf.py subsets the tandem repeat annotated VCF to only variants with tandem repeats identified  
    - filter_vcf_samps.py is used in filter_jasmine.nf to identify variants absent from the PanGenome and called by manta and delly  
* main.nf is the main nextflow script that is called to run the pipeline  
* docker_files/ contains the 5 dockerfiles needed to build the images used in the pipeline
* nextflow.config contains information for the different profiles; main configuration file for pipeline
* nextflow.local.config user-specific customization (?)
* nf_config.yaml config file with parameters defined 

### How to Run  
The nextflow pipeline expects as input a set of compressed, paired-end fastq files, the `nf_config.yaml` and `nextflow.local.config` files. It will output a single, concatenated VCF containing SNVs+INDELs called with Sentieon, CNVs called with Delly, SVs called with PanGenie, Manta, and Delly, and tandem repeat annotations output from Truvari. It also outputs a tarball containing the metrics files and the sorted BAM file with duplicates marked. The s3 bucket path is defined in `nf_config.yaml`

Below is a code block with an example for running the nextflow pipeline for WGS data. This is the default pipeline. 

```
nextflow /home/ec2-user/haffener_test_dir/sentieon_wgs/main.nf \
--proband_fastq1 /home/ec2-user/gencell_debug/V350204537_L03_13_433316_1.fq.gz \
--proband_fastq2 /home/ec2-user/gencell_debug/V350204537_L03_13_433316_2.fq.gz \
--platform ("CG", "MGI", "Element", or "Illumina") \
-c /home/ec2-user/haffener_test_dir/sentieon_wgs/nextflow.local.config \
-params-file /home/ec2-user/haffener_test_dir/sentieon_wgs/nf_config.yaml
```

The pipeline can also do secondary analysis for WES data by specifying the `--wes true` flag.  

```
nextflow /home/ec2-user/haffener_test_dir/sentieon_wgs/main.nf \
--proband_fastq1 /home/ec2-user/gencell_debug/V350204537_L03_13_433316_1.fq.gz \
--proband_fastq2 /home/ec2-user/gencell_debug/V350204537_L03_13_433316_2.fq.gz \
--platform ("CG", "MGI", "Element", or "Illumina") \
--target_bed /home/ec2-user/haffener_test_dir/GRCh38_protein_coding.bed \
-c /home/ec2-user/haffener_test_dir/sentieon_wgs/nextflow.local.config \
--wes true -params-file /home/ec2-user/haffener_test_dir/sentieon_wgs/nf_config.yaml
```

## Other Considerations
Regardless of the reference genome file specified in the params file, PanGenie calls against the GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set reference as it was used to generate the index files for PanGenie. 
For CNV calling with WES data, Sentieon states a limitation of the algorithm is that it was trained on diploid WGS samples. It is therefore recommended that with male samples, autosomes and sex chromosomes are called separately. See their example below:

```
 sentieon driver --interval AUTOSOMES_BED -t NUMBER_THREADS \
   -r REFERENCE -i DEDUPED_BAM --algo CNVscope \
   --model ML_MODEL/cnv.model TMP_DIPLOID_VCF

 sentieon driver --interval  -t NUMBER_THREADS \
   -r REFERENCE --algo CNVModelApply --model ML_MODEL/cnv.model \
   -v TMP_DIPLOID_VCF DIPLOID_VCF

sentieon driver --interval HAPLOID_BED -t NUMBER_THREADS \
   -r REFERENCE -i DEDUPED_BAM --algo CNVscope \
   --model ML_MODEL/cnv.model TMP_HAPLOID_VCF

 sentieon driver --interval  -t NUMBER_THREADS \
   -r REFERENCE --algo CNVModelApply --model ML_MODEL/cnv.model \
   -v TMP_HAPLOID_VCF HAPLOID_VCF
```