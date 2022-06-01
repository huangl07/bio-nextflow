#!/usr/bin/env nextflow
params.outdir = "demo"
params.help = false
params.method="GATK"
params.miss=0.3
params.maf=0.05
params.win=1000000
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
        file "pop.filtered.vcf.gz" into filter_vcf1,filter_vcf2,filter_vcf3,filter_vcf4,filter_vcf5,filter_vcf6
        file "pop.filtered.vcf.gz.tbi" into vcf_index
        file "chr.list" into chr_list,chr_list1
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
        file "group.list" into grouplist1,grouplist2
        file "xpclr.list" into xpclr_group,xpehh_group
        file "treemix.list" into treemix_group
    script:
    """
        perl ${baseDir}/bin/group.pl -i ${group} -o ./
        perl ${baseDir}/bin/remakegrolist.pl -in ${group} -out  treemix.list
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
        file "*.Tajima.D"  into tajimaD
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
        file "pixy_dxy.txt" into dxy
        file "pixy_fst.txt" into fst
        file "pixy_pi.txt" into pi
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
process merge_result{
    publishDir "${params.outdir}/05.merge_result", pattern:"*.detail"
    queue "DNA"
    cpus 1
    executor "slurm"
    memory "30G"
    cache 'lenient'
    input:
        file pi from pi
        file fst from fst
        file dxy from dxy
        file tajimaD from tajimaD.collect()
    output:
        path "*.detail" into merge_result_ch1, merge_result_ch2, merge_result_ch3,merge_result_ch4
    script:
    """
        Rscript ${baseDir}/bin//merge_SG_result.R --fst ${fst} --dxy ${dxy} --pi ${pi} --tajimaD ./ -m 8
    """
}

process fstthetapi{
    publishDir "${params.outdir}/06.fst-thetapi", pattern:"*"
    queue "DNA"
    cpus 1
    executor "slurm"
    memory "10G"
    cache 'lenient'
    input:
        file group from merge_result_ch1.flatten()
    output:
        file "*"
        file "*.select" into select_chr1
    script:
    """
        Rscript ${baseDir}/bin/fst_thetapi.R --merge ${group} --threshold 0.05

    """
}

process draw_fst_dxy{
    publishDir "${params.outdir}/07.fst-dxy", pattern:"*fst_dxy.pdf"
    queue "DNA"
    cpus 1
    executor "slurm"
    memory "10G"
    cache 'lenient'
    input:
        path group from merge_result_ch2.flatten()
    output:
        file "*"
        file "*.select" into select_ch2
    script:
    """
        Rscript ${baseDir}/bin/fst_dxy.R --merge ${group} --threshold 0.05
    """
}
process draw_fst_tajimaD_pi{
    publishDir "${params.outdir}/07.fst-dxy", pattern:"*fst_dxy.pdf"
    queue "DNA"
    cpus 1
    executor "slurm"
    memory "10G"
    cache 'lenient'
    input:
        path group from merge_result_ch3.flatten()
    output:
        path "*.pdf"
        path "*.png"
    script:
    """
        Rscript ${baseDir}/bin/fst-tajimaD-pi.R --merge ${group} --threshold 0.05
    """
}
//不知道为什么必须要这么写
xpclr_group.splitCsv(header:false,sep:'\t').combine(filter_vcf3).set{groups3}
//.map{row-> tuple(file(row[0]),file(row[1]),row(3))}.combine(filter_vcf3).into{groups1;groups2;groups3}
chr_list.splitCsv().combine(groups3).set{xpclr_groups}
//xpclr --input 01.filter/pop.filtered.vcf.gz  --samplesA 02.group/Cn.list --samplesB 02.group/DH2.list  --out ./Cn-DH2.xpclr --chr LG1 --size 100000
process xpclr{
   publishDir "${params.outdir}/08.xpclr", pattern:"*"
   queue "DNA"
   cpus 8
   executor "slurm"
   memory "30G"
   input:
        tuple chr,g1,g2,name,vcf from xpclr_groups
      // each chr from chrs
   output:
       file "*.xpclr" into  xpclr_reseult
   script:
   """
   xpclr --input ${vcf}  --samplesA ${g1} --samplesB ${g2}  --out ./${name}-${chr}.xpclr  --chr ${chr}  --size 100000
   """
}
process rearrange_xpclr_result{
    publishDir "${params.outdir}/09.xpclr", pattern:"*xpclr_abstract"
    queue "DNA"
    cpus 1
    executor "slurm"
    memory "30G"
    cache 'lenient'
    input:
        file xpclr from xpclr_reseult.collect()
    output:
        path "*xpclr.result" into xpclr_ch
    script:
    """
    ls *.xpclr|perl -ne '{chomp;@a=split("-",\$_);push @{\$stat{\$a[0]}},\$_;}END{foreach \$id(keys %stat){print "cat ",join(" ",@{\$stat{\$id}})," > \$id.xpclr.result\\n";}}'|sh
    """
}
process draw_manhattan_xpclr{
    publishDir "${params.outdir}/10.xpclrResult", pattern:"*"
    queue "DNA"
    cpus 1
    executor "slurm"
    memory "30G"
    cache 'lenient'
    input:
        path xpclr from xpclr_ch.flatten()
    output:
        file "*"
        file "*.select" into select_ch3
    script:
    """
         Rscript ${baseDir}/bin/xpclr.R --xpclr ${xpclr} --threshold 0.05
    """
}
//####################################################################################

grouplist2.splitCsv(header:false,sep:'\t').map{row-> tuple(row[0],file(row[1]))}.combine(filter_vcf6).set{group4}

chr_list1.splitCsv().combine(group4).set{chr_vcf_input}
process vcf_chr{
   publishDir "${params.outdir}/08.xpclr", pattern:"*"
   queue "DNA"
   cpus 8
   executor "slurm"
   memory "30G"
   input:
        tuple chr,gid,gtxt,vcf from chr_vcf_input
      // each chr from chrs
   output:
       file "${chr}-${gid}.vcf.gz" into  vcf_chr
   script:
   """
   bcftools view -r ${chr} -S ${gtxt} -m2 -M 2  --threads 8 $vcf > ${chr}-${gid}.vcf.gz
   """
}
process compare{
   publishDir "${params.outdir}/08.xpehh", pattern:"*"
   queue "DNA"
   cpus 8
   executor "slurm"
   memory "30G"
   input:
        file group from group_file
        file vcf from vcf_chr.collect()
   output:
        file "compare.list" into compare
   script:
   """
    perl ${baseDir}/bin/xpehh.pl -d ${group} -o compare.list
   """
}
compare.splitCsv(header:false,sep:'\t').map{row-> tuple(file(row[0]),file(row[1]),row[2],row[3])}.set{xpehh}

process selscan{
    publishDir "${params.outdir}/08.xpehh", pattern:"*"
   queue "DNA"
   cpus 8
   executor "slurm"
   memory "30G"
   input:
      tuple file(g1),file(g2),compare,chr  from xpehh
   output:
       file "${compare}-${chr}.xpehh.out" into  xpehh_result
   script:
   """
    plink --vcf ${g1} --recode --allow-extra-chr --double-id
    Rscript ${baseDir}/bin/map-xpehh.R --map plink.map --out total.map
    selscan --xpehh --vcf ${g1} --vcf-ref ${g2}  --map total.map --threads 8 --out ${compare}-${chr} --unphased --cutoff 0.05 --trunc-ok
    sed -i s/\\./${chr}/g ${compare}-${chr}.xpehh.out
   """
}


process rearrange_xpehh_result{
    publishDir "${params.outdir}/09.xpehh", pattern:"*xpehh_abstract"
    queue "DNA"
    cpus 1
    executor "slurm"
    memory "30G"
    cache 'lenient'
    input:
        file xpehh from xpehh_result.collect()
    output:
        path "*.xpehh.result" into xpehh_ch
    script:
    """
    ls *.xpehh.out|perl -ne '{chomp;@a=split("-",\$_);push @{\$stat{\$a[0]}},\$_;}END{foreach \$id(keys %stat){print "cat ",join(" ",@{\$stat{\$id}})," > \$id.xpehh.result\\n";}}'|sh
    """
}


process draw_manhattan_xpehh{
    publishDir "${params.outdir}/10.xpehhResult", pattern:"*"
    queue "DNA"
    cpus 1
    executor "slurm"
    memory "30G"
    cache 'lenient'
    input:
        path xpehh from xpehh_ch.flatten()
    output:
        file "*"
        file "*.select" into select_ch4
    script:
    """
         Rscript ${baseDir}/bin/xpehh.R --xpehh ${xpehh} --threshold 0.05
    """
}





















process prepairformix{
   publishDir "${params.outdir}/11.vcf2treemix", pattern:"*"
   queue "DNA"
   cpus 8
   executor "slurm"
   memory "30G"
   input:
       file vcf from filter_vcf4
       file list from treemix_group
      // each chr from chrs
   output:
       file "treemix/pop.tmix.gz" into treemix
       file "*"
   script:
   """
    mkdir treemix
    bgzip -d ${vcf}
    python ${baseDir}/bin/vcf2treemix.py -vcf pop.filtered.vcf -pop treemix.list -out treemix/
    gzip treemix/pop.tmix
   """
}
Channel.fromList([1,2,3,4,5]).combine(treemix).combine([1,2,3,4,5]).set{treemix1}

process treemix{
   publishDir "${params.outdir}/12.treemix", pattern:"*"
   queue "DNA"
   cpus 8
   executor "slurm"
   memory "30G"
   input:
      tuple m,file(treemix),n from treemix1
   output:
       file "pop.*" into mix_result
       file "*"
   script:
   """
    treemix -i ${treemix} -m ${m}  -o pop.${m}.${n} -k 1000 -bootstrap -global
   """
}

process optm{
   publishDir "${params.outdir}/13.treemix", pattern:"*"
   queue "DNA"
   cpus 8
   executor "slurm"
   memory "30G"
   input:
      file "*" from mix_result.collect()
   output:
       file "*"
   script:
   """
    Rscript  ${baseDir}/bin/treemix.R --input ./ --output treemix
   """
}


process RON{
   publishDir "${params.outdir}/14.RON", pattern:"*"
   queue "DNA"
   cpus 8
   executor "slurm"
   memory "30G"
   input:
      tuple vcf from filter_vcf5
      file group from group_file
   output:
       file "*"
   script:
   """
    plink  --vcf $vcf --make-bed --recode --out pop --double-id --allow-extra-chr
    perl ${baseDir}/bin/regroup.pl  -g ${group} -i pop.ped -o pop.group -m pop.map 
    Rscript ${baseDir}/bin/ROH.R --pop pop.group 
    rename "ROH - " ROH- *
   """
}