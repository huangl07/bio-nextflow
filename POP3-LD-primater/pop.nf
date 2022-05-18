#!/usr/bin/env nextflow
params.outdir = "demo"
params.help = false
params.method="GATK"
params.miss=0.3
params.maf=0.05
params.dep=5
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
    --dep   <num>   dep
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
        bcftools filter --threads 8 -i "F_Missing <=${params.miss} && MAF > ${params.maf} && FORMAT/DP < ${params.dep}" ${vcf}  > pop.filtered.vcf
    """
}
process vcf2tree{
    publishDir "${params.outdir}/02.vcf2tree", pattern:"*"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        file vcf from filter_vcf1
    output:
        file "*"
        file "pop.phylip" into phylip_file1,phylip_file2
        file "pop.fasta" into fasta_file1,fasta_file2
    script:
        
    """
        perl ${baseDir}/bin/vcf2tree.pl -i ${vcf} -o pop
    """
}

process modeltest{
    publishDir "${params.outdir}/03.modeltest", pattern:"*"
    queue "DNA"
    cpus 16
    executor "slurm"
    memory "150G"
    input:
        file phylip from phylip_file1
    output:
        file "pop.model.test.log" into model
        file "*"
    script:
        
    """
        modeltest-ng  -p 16 -d nt -i ${phylip}  -o pop.model.test
      
    """
}
process iqtree2{
    publishDir "${params.outdir}/03.modeltest", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "100G"
    input:
        file phylip from phylip_file2
        file model_file from model
    output:
        file "*"
    script:
    if(params.group){
        """
        less -S ${model_file}|grep iqtree|uniq|sed 's/>//g'|sed 's/iqtree/iqtree2 -st DNA -B 1000 -nt 8 -redo/g' > iqtree2.sh
        sh iqtree2.sh
        Rscript ${baseDir}/bin/tree.R --infile pop.phylip.treefile --outfile pop --group ${group_file}
        """
    }else{
        """
        less -S ${model_file}|grep iqtree|uniq|sed 's/>//g'|sed 's/iqtree/iqtree2 -st DNA -B 1000 -nt 8 -redo/g' > iqtree2.sh
        sh iqtree2.sh
        Rscript ${baseDir}/bin/tree.R --infile pop.phylip.treefile --outfile pop
        """
    }
}

process vcf2bed{
    publishDir "${params.outdir}/04.vcf2bed", pattern:"*"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        file vcf from filter_vcf2
    output:
        file "pop.*" into plink_files1,plink_file2
    script:
        
    """
        vcftools --vcf ${vcf} --plink --out pop
        plink --file pop --make-bed --out pop  
    """
}



process admixture{
    publishDir "${params.outdir}/05.admixture",pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        file plink from plink_files1.collect()
        each x from 1..20
    output:
        file "pop.*" into admixs
        file "*"
    script:
    """
        admixture pop.bed ${x} --cv -j8 > pop.${x}.log
        cut -f 1 pop.fam -d " "|paste - pop.${x}.Q > pop.${x}.result
    """
}
process CVresule{
    publishDir "${params.outdir}/06.CVresult",pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        file "*" from admixs.collect()
    output:
        file "*"
    script:
    """
       perl ${baseDir}/bin/admixture.pl -i ./ -o ./
    """
}
process pca{
    publishDir "${params.outdir}/06.pca", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "30G"
    input:
         file vcf from filter_vcf3
    output:
        file "pop.*" 
    script:

     if(params.group){
        """
        plink --vcf ${vcf} --pca --out pop --allow-extra-chr --double-id
        Rscript ${baseDir}/bin/pca.R --infile pop.eigenvec --outfile pop.pca --varfile pop.eigenval --group ${group_file}
        """
    }else{
        """
        plink --vcf ${vcf} --pca --out pop --allow-extra-chr --double-id
        Rscript ${baseDir}/bin/pca.R --infile pop.eigenvec --outfile pop.pca --varfile pop.eigenval
        """
    }
}
if(params.group){
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
    process LDdecay{
        publishDir "${params.outdir}/07.LDdecay", pattern:"*"
        queue "DNA"
        cpus 8
        executor "slurm"
        memory "30G"
        input:
            file vcf from filter_vcf4
            tuple gid,gfile from group1
        output:
            file "*"  
            file "${gid}.stat.gz" into ld_result
        script:

        """
        PopLDdecay --InVCF ${vcf} --OutStat ${gid} --MAF ${params.maf} --Miss ${params.miss} -SubPop ${gfile[0]}
        """
    }

    process LDdraw{
        publishDir "${params.outdir}/07.LDdecay", pattern:"*"
        queue "DNA"
        cpus 8
        executor "slurm"
        memory "30G"
        input:
            file $stat from ld_result.collect()
        output:
            file "*"  
        script:

        """
        ls *.stat.gz|perl -ne 'chomp;\$a=\$_;\$a=~s/\\.stat\\.gz//g;print \$a,"\\t",\$_,"\\n"' > pop.list
        Rscript ${baseDir}/bin/ld-decay.R --list pop.list --outfile pop
        """
    }
}else{

    process LDdecaynogroup{
        publishDir "${params.outdir}/07.LDdecay", pattern:"*"
        queue "DNA"
        cpus 8
        executor "slurm"
        memory "30G"
        input:
            file vcf from filter_vcf4
        output:
            file "*"  
        script:
            """
            PopLDdecay --InVCF ${vcf} --OutStat pop --MAF 0.05 --Miss 0.3
            Rscript ${baseDir}/bin/ld-decay.R --infile pop.stat.gz --outfile pop.ld
            """
    }
}


    process primater{
    publishDir "${params.outdir}/08.parameter", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "100G"
    input:
         file vcf from filter_vcf5
         file plink from plink_file2.collect()

    output:
        file "*" 
    script:
       if(params.group){
        """
        populations -V ${vcf} -M ${group_file} --fstats -O ./ -t 8
        """
        }else{
        """
        cut -f 1 pop.fam -d " "|perl -ne 'chomp;print \$_,"\\tg1\\n";' > group.list
        populations -V ${vcf} -M group.list -O ./ -t 8
        """
    }       
}


workflow.onComplete {
    println "Pipeline completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
