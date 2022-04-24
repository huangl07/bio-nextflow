#!/usr/bin/env nextflow
params.outdir = "demo"
params.winsize="1000000"
params.stepsize="10000"
params.help = false
params.popt="F2"
params.bulksize=30
params.bootstrap=1000
params.pvalue=0.05
params.grade=false
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
   

    --vcf   <file>  input vcf file
    --outdir   <dir>   output dir
    --chr   <file>  chr list for draw
    --group <file>   group file
        mbid must be given
        
        wpid    flag    phdep   pldep
        mpid    flag    phdep   pldep

    --popt  <str>   population 
    ###############################################################
    --bulksize  <num>   bulk size
    --winsize <num>     windows size
    --stepsize <num>    step size
    ###############################################################
    --bootstrap <num>   bootstrap number
    --pvalue    <num>   pvalue
    ################################################################
    --grade <bool>  do grade for multi bulk than 2 

    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}
vcf=file(params.vcf)
group=file(params.group)
chr=file(params.chr)

if(!params.grade){
process vcf2table{
    publishDir "${params.outdir}/01.vcf2table", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file vcf from vcf
        file group from group
    output:
        file "pop.index" into index1,index2,index3
        file "*"
    script:
        """
        perl ${baseDir}/bin/vcf2table.pl --vcf ${vcf} --out pop.table --group ${group} -popt ${params.popt} -vtype ALL
        perl ${baseDir}/bin/bsa_calc.pl --table pop.table --group ${group} --popt ${params.popt} -out pop.index
        """
}
process indexslid{
    publishDir "${params.outdir}/02.index-slid", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file index_file from index1
    output:
        file "*"
    script:
    """

        Rscript ${baseDir}/bin/slidingwin-index.R --infile ${index_file} --outfile pop --winsize ${params.winsize} --stepsize  ${params.stepsize}  --method bp
        Rscript ${baseDir}/bin/bootstrap-index.R --infile pop.sliding.detail --bulksize ${params.bulksize} --outfile pop --bootstrap ${params.bootstrap} --popstruc ${params.popt} --qvalue ${params.pvalue}
        Rscript ${baseDir}/bin/manhattan-index.R --result pop.bootstrap.result --chr ${chr} --output pop.ED  --pcol slidingED --threshold 0.9995 --xlab chromosome --ylab "Eudiance Metric" 
        Rscript ${baseDir}/bin/manhattan-index.R --result pop.bootstrap.result --chr ${chr} --output pop.index --pcol delta --lcol slidingD --qcol CI --xlab chromosome --ylab "delta-index"
        Rscript ${baseDir}/bin/region.R --infile pop.bootstrap.result --ccol X.chr --pos pos --loess slidingD --CI CI --outfile pop.index.region --number 10 
        Rscript ${baseDir}/bin/region.R --infile pop.bootstrap.result --ccol X.chr --pos pos --loess slidingED --quantile 0.995 --outfile pop.ED.region --number 10 

    """
}
process loesscalc{
    publishDir "${params.outdir}/04.loess", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file table_file from index2
    output:
        file "*"
    script:
    """
        
        Rscript ${baseDir}/bin/loess-index.R --infile ${table_file} --out pop.index.loess
        Rscript ${baseDir}/bin/bootstrap-index.R --infile pop.index.loess --bulksize ${params.bulksize} --outfile pop --bootstrap ${params.bootstrap} --popstruc ${params.popt} --qvalue ${params.pvalue}
        Rscript ${baseDir}/bin/manhattan-index.R --result pop.bootstrap.result --pcol delta --qcol CI --lcol loess --xlab chromosome --ylab "loess" --output loess --chr ${chr}
        Rscript ${baseDir}/bin/region.R --infile pop.bootstrap.result --ccol X.chr --pos pos --loess loess --CI CI --outfile pop.loess.region --number 10 

    """
}
process Gcalc{
    publishDir "${params.outdir}/03.Gprime", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file index_file from index3
    output:
        file "*"
    script:
    """
        Rscript ${baseDir}/bin/G-calc.R --infile ${index_file} --outfile pop.Gprime --winsize ${params.winsize} --method bp
        Rscript ~dna/Pipeline/NextFlow/bio-nextflow/BSA/bin/manhattan-index.R --result pop.Gprime --pcol G --lcol Gprime --xlab chromosome --ylab "Gprime" --output Gprime  --chr ${chr}
        Rscript ${baseDir}/bin/region.R --infile pop.Gprime --ccol X.chr --pos pos --loess Gpval --CI threshold --outfile pop.Gprime.region --number 10 

    """
}
}else{
process params.grade{
    publishDir "${params.outdir}/01.vcf2table", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file vcf from vcf
    output:
        file "*"
    script:
    """
        perl ${baseDir}/bin/vcf2table_multi.pl --vcf ${vcf} --out pop.table --group ${group} -popt ${params.popt} -vtype ALL
        Rscript ${baseDir}/bin/ridit.R --index pop.table pop.table --out pop --window ${params.winsize} --pvalue ${params.pvalue} --group ${params.group}
        Rscript ${baseDir}/bin/manhattan-index.R --result pop.denoise.result --pcol logP --qcol 1 --xlab chromosome --ylab "loess" --output loess --chr ${chr}
        Rscript ${baseDir}/bin/region.R --infile pop.denoise.result --ccol X.chr --pos pos --loess loess --CI 1 --outfile pop.loess.region --number 10 
    """
}
}
