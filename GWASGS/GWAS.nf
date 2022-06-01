#!/usr/bin/env nextflow
params.outdir = "demo"
params.help = false
params.miss=0.3
params.maf=0.05
params.dep=5
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --vcf     <file>   output dir
    --trt     <file>  input group file    
    --outdir  <dir>   output dir

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

vcf=file(params.vcf)
trt =file(params.trt)


process vcffilter{
    publishDir "${params.outdir}/01.filter", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
        cache 'lenient'

    input:
        file vcf from vcf
    output:
        file "pop.filtered.vcf" into filter_vcf1,filter_vcf2
        file "pop.hmp.txt" into hmp_files
    script:
    """
        bcftools filter --threads 8 -O v -i "F_Missing <=${params.miss} && MAF > ${params.maf} && FORMAT/DP < ${params.dep}" ${vcf}  > pop.filtered.vcf
        perl ${baseDir}/bin/vcf2hapmap.pl -i pop.filtered.vcf -o  pop.hmp.txt
    """
}

process blocks{
    publishDir "${params.outdir}/02.blocks", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    input:
        file vcf from filter_vcf2
    output:
        file "pop.blocks" into blocks
    script:
    """
        plink --vcf ${vcf} --out pop  --blocks no-pheno-req --allow-extra-chr --double-id
    """
}
process trtprepair{
    publishDir "${params.outdir}/02.trt", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
        cache 'lenient'

    input:
        file trt from trt
    output:
        file "trt.list" into trt_list1,trt_list2
    script:
    """
        perl ${baseDir}/bin/trt.pl -i ${trt} -o ./
    """
}
trt_list1.splitCsv(header:false,sep:'\t').combine(filter_vcf1).into{gwas}
trt_list2.splitCsv(header:false,sep:'\t').combine(hmp_files).into{gs}
process GWAS{
    publishDir "${params.outdir}/03.GWAS", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    input:
        tuple trtid,trt,vcf from gwas
    output:
        tuple trtid,"GWAS.result.csv" into gwasresult
        file "*signals.csv" into signals
        file "*"
    script:
    """
        Rscript ${baseDir}/bin/MVP.single.R --vcf ${vcf} --trait ${trt} --output ./ --pvalue 0.01
    """
}

process GS{
    publishDir "${params.outdir}/03.GWAS", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "150G"
    input:
        tuple trtid,trt,vcf from gs
    output:
        file "*"
    script:
    """
        Rscript ${baseDir}/bin/rrblup.R --vcf ${vcf} --trait ${trt} --output ./ 
    """
}