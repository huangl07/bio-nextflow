#!/usr/bin/env nextflow
params.outdir = "demo"
params.help = false
params.method="GATK"
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --vcf       <file>  input vcf file
    --outdir    <dir>   output dir
    --group     <file>  input group file    
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
        file "pop.filtered.vcf" into filter_vcf1,filter_vcf2,filter_vcf3
    script:
        
    """
        bcftools filter --threads 8  -i "F_Missing <=0.2 && MAF > 0.05" ${vcf}  > pop.filtered.vcf
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
    cpus 8
    executor "slurm"
    memory "100G"
    input:
        file phylip from phylip_file1
    output:
        file "pop.model.test.log" into model
        file "*"
    script:
        
    """
        modeltest-ng  -p 8 -d nt -i ${phylip}  -o pop.model.test
      
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
    if(params.group)
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


num=Channel.from(2..20)

process admixture{
    publishDir "${params.outdir}/05.admixture",pattern:"*"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        file plink from plink_files1.collect()
        val x from num
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
        file "pop.*" into plink_files
    script:

     if(params.group)
        """
        plink --vcf ${vcf} --pca --out pop --allow-extra-chr 
        Rscript ${baseDir}/bin/pca.R --infile pop.eigenvec --outfile pop.pca --varfile pop.eigenval --group ${group_file}
        """
    }else{
        """
        plink --vcf ${vcf} --pca --out pop --allow-extra-chr 
        Rscript ${baseDir}/bin/pca.R --infile pop.eigenvec --outfile pop.pca --varfile pop.eigenval
        """
    }
        
    
}