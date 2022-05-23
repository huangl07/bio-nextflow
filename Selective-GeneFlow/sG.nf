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
        cache 'lenient'

    input:
        file vcf from vcf_file
    output:
        file "pop.filtered.vcf.gz" into filter_vcf1,filter_vcf2,filter_vcf3
        file "pop.filtered.vcf.gz.tbi" into vcf_index
        file "chr.list" into chr_list
    script:
    """
        bcftools filter --threads 8  -O z -i "F_Missing <=${params.miss} && MAF > ${params.maf} && FORMAT/DP < ${params.dep}" ${vcf}  > pop.filtered.vcf.gz
        tabix pop.filtered.vcf.gz
        less -S pop.filtered.vcf.gz|grep "#"|perl -ne 'if(/##contig=<ID=([^,]*)/){print \$1,"\\n"}' > chr.list
    """
}


process group_preparie{
    publishDir "${params.outdir}/02.group", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
        cache 'lenient'

    input:
        file group from group_file
    output:
        file "group.list" into grouplist1
        file "xpclr.list" into grouplist2
    script:
    """
        perl ${baseDir}/bin/group.pl -i ${group} -o ./
    """
}
grouplist1.splitCsv(header:false,sep:'\t').map{row-> tuple(row[0],file(row[1]))}.combine(filter_vcf1).set{groups}

process tajimaD1{
    publishDir "${params.outdir}/03.tajimaD", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    cache 'lenient'
    input:
        tuple gid,file(gfile),file(vcf) from groups
    output:
        file "*.Tajima.D" 
    script:
    """
       vcftools --gzvcf ${vcf} --remove-indels --TajimaD ${params.win} --keep ${gfile} --out ${gid}
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
        file index from vcf_index
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

grouplist2.splitCsv(header:false,sep:'\t').combine(filter_vcf3).set{groups3}
//.map{row-> tuple(file(row[0]),file(row[1]),row(3))}.combine(filter_vcf3).into{groups1;groups2;groups3}
chr_list.splitCsv().combine(groups3).set{xpclr_group}
//xpclr --input 01.filter/pop.filtered.vcf.gz  --samplesA 02.group/Cn.list --samplesB 02.group/DH2.list  --out ./Cn-DH2.xpclr --chr LG1 --size 100000
process xpclr{
   publishDir "${params.outdir}/05.xpclr", pattern:"*"
   queue "DNA"
   cpus 8
   executor "slurm"
   memory "30G"
   input:
       tuple chr,g1,g2,name,vcf from xpclr_group
      // each chr from chrs
   output:
       file "*"
   script:
   """
   xpclr --input ${vcf}  --samplesA ${g1} --samplesB ${g2}  --out ./${name}-${chr}.xpclr  --chr ${chr}  --size 100000
   """
}