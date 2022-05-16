#!/usr/bin/env nextflow
params.outdir = "demo"
params.help = false
params.method="GATK"
params.miss=0.3
params.maf=0.05
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --vcf       <file>  input vcf file
    --outdir    <dir>   output dir
    --group     <file>  input group file    
    --miss  <num>   miss
    --maf   <num>   maf
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

vcf_file = files(params.vcf)
if(params.group){
    group_file = file(params.group)
}
process vcffilter{
    publishDir "${params.outdir}/01.filter", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    input:
        file vcf from vcf_file
    output:
        file "*"
        file "pop.filtered.vcf" into filter_vcf1,filter_vcf2,filter_vcf3,filter_vcf4,filter_vcf5
    script:
        
    """
        bcftools filter --threads 8  -i "F_Missing <=${params.miss} && MAF > ${params.maf}" ${vcf}  > pop.filtered.vcf
    """
}