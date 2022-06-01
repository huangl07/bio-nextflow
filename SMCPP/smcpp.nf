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

    --vcf   <file>  normal bam file
    --group    <file>  tumor  bam file
    --outdir   <dir>   output dir
    --miss  <num>   miss
    --maf   <num>   maf
    --dep   <num>   dep
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
        file "pop.vcf.gz" into filter_vcf1,filter_vcf2,filter_vcf3
        file "pop.vcf.gz.tbi" into vcf_index
        file "chr.list" into chr_list
    script:
    """
        bcftools filter --threads 8  -O z -i "F_Missing <=${params.miss} && MAF > ${params.maf} && FORMAT/DP < ${params.dep}" ${vcf}  > pop.filtered.vcf.gz
        plink --vcf pop.filtered.vcf.gz --keep-allele-order --recode vcf-iid --out pop --threads 8 --double-id --allow-extra-chr
        bgzip pop.vcf
        tabix pop.vcf.gz
        less -S pop.vcf.gz|grep "#"|perl -ne 'if(/##contig=<ID=([^,]*)/){print \$1,"\\n"}' > chr.list
    """
}

process vcf2psmc{
    publishDir "${params.outdir}/02.vcf2psmc", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    input:
        file group from group_file
        file vcf from filter_vcf1
    output:
        file ("psmc.list") into psmclist
        file "*"
    script:
    """
        perl ${baseDir}/bin/vcf2psmc.pl -vcf ${vcf} -out ./ -gro ${group}
    """
}

psmclist.splitCsv(header:false,sep:'\t').set{psmcs}

process psmc{
    publishDir "${params.outdir}/03.psmc", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    input:
        tuple popid,sampleID,psmc from psmcs
    output:
        file "${popid}-${sampleID}.psmc" into psmc_result
    script:
    """
        psmc -N25 -t15 -r5 -p \"4+25*2+4+6\" -o ${popid}-${sampleID}.psmc ${psmc}
    """
}

process drawpsmc{
    publishDir "${params.outdir}/04.drawpsmc", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    input: 
        file psmc from psmc_result.collect()
    output:
        file "*"
    script:
    """
    ls *.psmc|perl -ne '{chomp;@a=split("-",\$_);push @{\$stat{\$a[0]}},\$_}END{foreach \$id(keys %stat){print "cat ",join(" ",@{\$stat{\$id}}),"> \$id.psmc.result\\n"}}' |sh
    ls *.psmc.result|perl -ne 'chomp;@a=split(/\\./,\$_);print \$a[0],"\\t",\$_,"\\n"' > result.list
    Rscript ${baseDir}/bin/drawPSMC.R --infile result.list --outdir ./
    """
}





process vcf2hamp{
    publishDir "${params.outdir}/07.vcf2hamp", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
        cache 'lenient'
    input:
        file vcf from filter_vcf2
    output:
        file "pop.msmc.txt" into msmc
    script:
    """
        perl  ${baseDir}/bin/vcf2hamp-ms.pl  -i ${vcf} -o pop.msmc.txt 
    """
}
//
process group{
    publishDir "${params.outdir}/08.group2msmc", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    cache 'lenient'
    input:
        file group from group_file
    output:
        file "group.list" into group_list1,group_list2,group_list3
    script:
    """
        perl ${baseDir}/bin/group.pl -group ${group} -out group.list
    """
}

group_list1.splitCsv(header:false,sep:'\t').combine(msmc).set{msmc}
process hamp2msmc{
    publishDir "${params.outdir}/09.hamp2msmc", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "500G"
    cache 'lenient'
    input:
        tuple pop,sampleID,file(msmc)  from msmc 
    output:
        file "${pop}.msmc.result" into msmc_result
    script:
    """
        python ${baseDir}/bin/make_input_MSMC_from_callsTab.py -i ${msmc} -o ${pop}.msmc -s ${sampleID} -m 1
        msmc2 --fixedRecombination --nrThreads=8 -o ${pop}.msmc.result *_${pop}.msmc
    """
}
process drawmsmc{
    publishDir "${params.outdir}/10.drawmsmc", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    input: 
        file msmc from msmc_result.collect()
    output:
        file "*"
    script:
    """
    ls *.msmc.result|perl -ne 'chomp;@a=split(/\\./,\$_);print \$a[0],"\\t",\$_,"\\n"' > result.list
    Rscript ${baseDir}/bin/drawmsmc.R --infile result.list --outdir ./
    """
}



group_list2.splitCsv(header:false,sep:'\t')
              .combine(filter_vcf3).combine(vcf_index).set{smc_group}
chr_list.splitCsv(header:false,sep:'\t').combine(smc_group).set{smcvcf}
process vcf2smc{
    publishDir "${params.outdir}/11.smcpp", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    cache 'lenient'
    input:
        tuple chr,pop,sid,file(vcf),file(index) from smcvcf
    output:
        file("${pop}-${chr}.smc.gz") into smc
        file "*"
    script:
    """
       smc++ vcf2smc --cores 8 --ignore-missing --missing-cutoff 100 ${vcf} ${pop}-${chr}.smc.gz ${chr} ${pop}:${sid}
    """
}
//
group_list3.splitCsv(header:false,sep:'\t').set{groups}
process smcpp{
    publishDir "${params.outdir}/12.smcpp", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    //clusterOptions '--nodelist=compute-3-112,compute-3-113' 
    container '10.2.4.236:5000/smcpp'
    input:
        file smc from smc.collect()
        tuple pop,gid from groups
    output:
        file "${pop}.model.json" into json
         file "*"
    script:
    """
        mkdir ${pop}
        cp `readlink -f ${pop}-*.smc.gz` ${pop}
        cd ${pop}
        docker run --rm -v \$PWD:/mnt 10.2.4.236:5000/smcpp cv --cores 8 -o ./ 1.25e-8 ${pop}-*.smc.gz
        mv model.final.json ../${pop}.model.json
    """
}

process drawsmc{
    publishDir "${params.outdir}/13.smcpp", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
        cache 'lenient'
    input:
        file json from json.collect()
    output:
        file "*"
    script:
    """
        smc++ --cores 8  plot pop.pdf *.json
    """
}