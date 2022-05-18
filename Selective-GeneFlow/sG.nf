#!/usr/bin/env nextflow
params.outdir = "demo"
params.help = false
params.method="GATK"
params.miss=0.3
params.maf=0.05
params.win=10000
params.dep=7
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
    --win   <num>   window size
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

vcf_file = files(params.vcf)
group_file = files(params.group)

process vcffilter{
    publishDir "${params.outdir}/01.filter", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    input:
        file vcf from vcf_file
    output:
        file "pop.filtered.vcf.gz" into filter_vcf1,filter_vcf2
        file "pop.filtered.vcf.gz.tbi" into index1,index2
        file "*"
    script:
        
    """
        bcftools filter --threads 8  -O z -i "F_Missing <=${params.miss} && MAF > ${params.maf} && FORMAT/DP < ${params.dep}" ${vcf}  > pop.filtered.vcf.gz
        tabix pop.filtered.vcf.gz
    """
}


process group_preparie{
    publishDir "${params.outdir}/02.group", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    input:
        file group from group_file
    output:
        file "group.list" into grouplist1,grouplist2
        file "*"
    script:
    """
        perl ${baseDir}/bin/group.pl -i ${group} -o ./
    """
}

grouplist1.splitCsv(header:false,sep:'\t').groupTuple().set{group1}
grouplist2.splitCsv(header:false,sep:'\t').groupTuple().set{group2}


process tajimaD{
    publishDir "${params.outdir}/03.tajimaD", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    input:
        file vcf from filter_vcf1
        tuple gid,gfile from group2
        file index from index1
    output:
        file "*"
    script:
    """
       vcftools --gzvcf ${vcf} --remove-indels --keep ${gfile[0]} --out ${gid} --TajimaD ${params.win}
    """
}

process pixy{
    publishDir "${params.outdir}/04.pixD", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    input:
        file vcf from filter_vcf2
        file index from index2
        file group from group_file
    output:
        file "*"
    script:
    """
        pixy --stats pi fst dxy \
        --vcf ${vcf} \
        --populations ${group}\
        --window_size 10000 \
        --n_cores 8 \
        --bypass_invariant_check yes
    """
}