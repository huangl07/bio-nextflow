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
        file "pop.filtered.vcf" into filter_vcf1,filter_vcf2,filter_vcf3,filter_vcf4,filter_vcf5,filter_vcf6,filter_vcf7,filter_vcf8,filter_vcf9
    script:
        
    """
        bcftools filter --threads 8 -i "F_Missing <=${params.miss} && MAF > ${params.maf} && MIN(FORMAT/DP) > ${params.dep}" ${vcf}  > pop.filtered.vcf
    """
}
if(params.group){
    process group{
        publishDir "${params.outdir}/12.group-split", pattern:"*"
        queue "DNA"
        cpus 8
        executor "slurm"
        memory "30G"
        input:
            file group from group_file
        output:
            file "*"
            file "id.list" into group_split
        script:
        """
            less -S ${group}|perl -ne '{@a=split;push @{\$stat{\$a[1]}},\$a[0]}END{foreach \$id(keys %stat){open Out,">\$id.list";print Out join("\n",@{\$stat{\$id}});close Out;print \$id,"\t",`readlink -f \$id.list`}}'  > id.list
        """
    }
    group_split.splitCsv(header:false,sep:"\t").combine(filter_vcf9).set{nneinput}
    process effpopulation1{
        publishDir "${params.outdir}/13.Ne", pattern:"*"
        queue "DNA"
        cpus 8
        executor "slurm"
        memory "30G"
        input:
            tuple gid,group,vcf from nneinput
        output:
            file "*"
        script:
        """ 
            bcftools view  -S ${group} -O v  ${vcf} --threads 8 > ${gid}.vcf
            plink --vcf ${gid}.vcf --recode --out ${gid} --double-id  --allow-extra-chr 
            SNeP1.1 -threads 8 -ped ${gid}.ped
        """
    }
}else{
        process effpopulation2{
        publishDir "${params.outdir}/13.Ne", pattern:"*"
        queue "DNA"
        cpus 8
        executor "slurm"
        memory "30G"
        input:
            file vcf from filter_vcf9
        output:
            file "*"
        script:
        """
            plink --vcf ${vcf} --recode --out pop --double-id  --allow-extra-chr
            SNeP1.1 -threads 8 -ped pop
        """
        }
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


    process parameter{
    publishDir "${params.outdir}/08.parameter", pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "100G"
    input:
         file vcf from filter_vcf5
            file group from group_file

    output:
        file "*" 
    script:
       if(params.group){
        """
        populations -V ${vcf} -M ${group} --fstats -O ./ -t 8
        less -S  pop.filtered.p.sumstats_summary.tsv|perl -ne 'chomp;s/#//g;last if(/variant and fixed/);print \$_,"\n" if(!/Variant/);' > stack.list
        perl ${baseDir}/bin/pic.pl -i ${vcf} -g ${group} -o pic
        Rscript ${baseDir}/bin/diversity.R --vcf ${vcf} --group ${group} --out ./ --pic pic.stat --stack stack.list
        """
        }else{
        """
        cut -f 1 pop.fam -d " "|perl -ne 'chomp;print \$_,"\\tg1\\n";' > group.list
        populations -V ${vcf} -M group.list -O ./ -t 8
        """
        }       
    }
    process Amova{
        publishDir "${params.outdir}/09.Amova", pattern:"*"
        queue "DNA"
        cpus 8
        executor "slurm"
        memory "100G"
        input:
            file vcf from filter_vcf6
            file group from group_file

        output:
            file "*" 
        script:
            if(params.group){
                """
                Rscript ${baseDir}/bin/Amova.R --vcf ${vcf} --group ${group} -o pop 
                """
            }
    }
 
    process ibc{
        publishDir "${params.outdir}/10.ibc", pattern:"*"
        queue "DNA"
        cpus 8
        executor "slurm"
        memory "100G"
        input:
            file vcf from filter_vcf7
        output:
            file "*" 
        script:
        """
        plink --vcf ${vcf} --het --ibs-matrix --double-id --out pop --allow-extra-chr
        """    
    }
    process ibd{
        publishDir "${params.outdir}/11.idb", pattern:"*"
        queue "DNA"
        cpus 8
        executor "slurm"
        memory "100G"
        input:
            file vcf from filter_vcf8
        output:
            file "*" 
        script:
        """
        java -jar $prefix/java-jar/beagle.27Jan18.7e1.jar gt=${vcf}  out=pop ibd=true impute=false
        """    
    }

workflow.onComplete {
    println "Pipeline completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
