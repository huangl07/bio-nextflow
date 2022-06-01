#!/usr/bin/env nextflow
params.outdir = "demo"
params.winsize="1000000"
params.stepsize="10000"
params.help = false
params.popt="F2"
params.bulksize=30
params.bootstrap=1000
params.pvalue=0.001
params.grade=false
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
   

    --vcf   <file>  input vcf file
    --outdir   <dir>   output dir
    --chr   <file>  chr list for draw
    ##############################################################
    --gff   <file>  gff file
    --GOanno    <file>  GOanno
    --KEGGanno  <file>  KEGGanno
    ##############################################################
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
    --mutmap

    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}
vcf=file(params.vcf)
group=file(params.group)
chr=file(params.chr)
gff_file=file(params.gff)
GOanno=file(params.GOanno)
KEGGanno=file(params.KEGGanno)
if(!params.grade){
process vcf2table{
    publishDir "${params.outdir}/01.vcf2table", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file vcf from vcf
        file group from group
    output:
        file "pop.index" into index1,index2,index3,index4
        file "pop.table" into table_file
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
        tuple val("index"),file("pop.index.region") into index_region
        file "*"
    script:
    if(params.mutmap){
    """

        Rscript ${baseDir}/bin/slidingwin-index.R --infile ${index_file} --outfile pop --winsize ${params.winsize} --stepsize  ${params.stepsize}  --method bp
        Rscript ${baseDir}/bin/bootstrap-index.R --infile pop.sliding.detail --bulksize ${params.bulksize} --outfile pop --bootstrap ${params.bootstrap} --popstruc ${params.popt} --qvalue ${params.pvalue} 
        Rscript ${baseDir}/bin/manhattan-index.R --result pop.bootstrap.result --chr ${chr} --output pop.index --pcol delta --lcol slidingD --qcol CI --xlab chromosome --ylab "delta-index"
        Rscript ${baseDir}/bin/region.R --infile pop.bootstrap.result --ccol X.chr --pos pos --loess slidingD --CI CI --outfile pop.index.region --number 10 
    """
    }else{
        """
        Rscript ${baseDir}/bin/slidingwin-index.R --infile ${index_file} --outfile pop --winsize ${params.winsize} --stepsize  ${params.stepsize}  --method bp
        Rscript ${baseDir}/bin/bootstrap-index.R --infile pop.sliding.detail --bulksize ${params.bulksize} --outfile pop --bootstrap ${params.bootstrap} --popstruc ${params.popt} --qvalue ${params.pvalue}
        Rscript ${baseDir}/bin/manhattan-index.R --result pop.bootstrap.result --chr ${chr} --output pop.index --pcol delta --lcol slidingD --qcol CI --xlab chromosome --ylab "delta-index"
        Rscript ${baseDir}/bin/region.R --infile pop.bootstrap.result --ccol X.chr --pos pos --loess slidingD --CI CI --outfile pop.index.region --number 10 
        """
    }
}
process loesscalc{
    publishDir "${params.outdir}/04.loess", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file table_file from index2
    output:
        file "*"
        tuple val("loess"),file("pop.loess.region") into loess_region

    script:
    """
        
        Rscript ${baseDir}/bin/loess-index.R --infile ${table_file} --out pop.index.loess
        Rscript ${baseDir}/bin/bootstrap-index.R --infile pop.index.loess --bulksize ${params.bulksize} --outfile pop --bootstrap ${params.bootstrap} --popstruc ${params.popt} --qvalue ${params.pvalue}
        Rscript ${baseDir}/bin/manhattan-index.R --result pop.bootstrap.result --pcol delta --qcol CI --lcol loess --xlab chromosome --ylab "loess" --output loess --chr ${chr}
        Rscript ${baseDir}/bin/region.R --infile pop.bootstrap.result --ccol X.chr --pos pos --loess loess --CI CI --outfile pop.loess.region --number 10 

    """
}
    if(!params.mutmap){
        process Gcalc{
            publishDir "${params.outdir}/03.Gprime", pattern:"*"
            queue "DNA"
            executor "slurm"
            input:
                file index_file from index3
            output:
                file "*"
                tuple val("Gprime"),file("pop.Gprime.region") into Gprime_region
            script:
            """
                Rscript ${baseDir}/bin/G-calc.R --infile ${index_file} --outfile pop.Gprime --winsize ${params.winsize} --method bp
                Rscript ~dna/Pipeline/NextFlow/bio-nextflow/BSA/bin/manhattan-index.R --result pop.Gprime --pcol G --lcol Gprime --qcol GprimeT --xlab chromosome --ylab "Gprime" --output Gprime  --chr ${chr}
                Rscript ${baseDir}/bin/region.R --infile pop.Gprime --ccol X.chr --pos pos --loess Gprime --CI GprimeT --outfile pop.Gprime.region --number 10 

            """
        }
        process EDslid{
        publishDir "${params.outdir}/02.index-slid", pattern:"*"
        queue "DNA"
        executor "slurm"
        input:
            file index_file from index1
        output:
            file "*"
            tuple val("ED"),file("pop.ED.region") into ED_region
        script:
        """
        Rscript ${baseDir}/bin/slidingwin-index.R --infile ${index_file} --outfile pop --winsize ${params.winsize} --stepsize  ${params.stepsize}  --method bp
        Rscript ${baseDir}/bin/manhattan-index.R --result pop.sliding.result --chr ${chr} --output pop.ED  --pcol slidingED --threshold 0.9995 --xlab chromosome --ylab "delta-index"
        Rscript ${baseDir}/bin/region.R --infile pop.sliding.result --ccol X.chr --pos pos --loess slidingED --quantile 0.995 --outfile pop.ED.region --number 10 
        """
        }
    }


    Channel.from(index_region).mix(loess_region).set{regions}
    if(!params.mutmap){
        Channel.from(regions).mix(ED_region,Gprime_region).set{region}
    }


}else{
process params.grade{
    publishDir "${params.outdir}/01.vcf2table", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file vcf from vcf
        file group from group
    output:
        file "*"
        tuple val("ridit"),"pop.ridit.region" into region
        file "pop.table" into table_file
    script:
    """
        perl ${baseDir}/bin/vcf2table_multi.pl --vcf ${vcf} --out pop.table --group ${group} -popt ${params.popt} -vtype ALL
        Rscript ${baseDir}/bin/ridit.R --index pop.table --out pop --window ${params.winsize} --pvalue ${params.pvalue} --group ${params.group}
        Rscript ${baseDir}/bin/manhattan-index.R --result pop.denoise.result --pcol logP --threshold 0.999 --xlab chromosome --ylab "loess" --output loess --chr ${chr}
        Rscript ${baseDir}/bin/region.R --infile pop.denoise.result --ccol X.chr --pos pos --loess loess --quantile 0.999 --outfile pop.ridit.region --number 10 
    """
}
}

process Orgdb {
    publishDir "${params.outdir}/00.orgdb", pattern:"*"
    executor 'slurm'
    queue 'DNA'
    cpus 2
    cache 'lenient'
    input:
        file GO_anno from GOanno
        file KEGG_anno from KEGGanno
    output:
        tuple path("TERM2GENE.txt"),path("TERM2NAME.txt") into term_ch
        path(LIB) into EnrichDb
    script:
        """
        mkdir LIB
        awk -v RS='----+' 'NR==2' ${KEGGanno} | awk -v RS="////" '{print}' | awk NF > KEGG_raw.result
        python  ${baseDir}/bin//KEGG_rearrange.py -filename KEGG_raw.result -outfile KEGG_annotation.csv
        cat KEGG_annotation.csv | awk -F "\\t" 'NR>1{print \$4"\\t"\$1}' | grep -v "^-" > TERM2GENE.txt
        cat KEGG_annotation.csv | awk -F "\\t" 'NR>1{print \$4"\\t"\$3}' | grep -v "^-" > TERM2NAME.txt
        Rscript ${baseDir}/bin//make_orgdb.R --gofile ${GOanno} --keggfile KEGG_annotation.csv --outpath .
        R CMD INSTALL org.Ddemo.eg.db --library=./LIB/
        """
} 
process params.enrich{
    publishDir "${params.outdir}/05.enrich", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        tuple method,region from region
        tuple path(term2gene),path(term2name) from term_ch
        file gff from gff_file
        path LIB from EnrichDb
        file GO_anno from GOanno
        file table from table_file
        file KEGG_file from KEGGanno
    output:
        file "*"
    script:
    if(region.countLines() != 0){
    """
        mkdir ${method}
        cd ${method}
        grep mRNA ../${gff} | awk '{print \$1"\\t"\$4"\\t"\$5"\\t"\$9}' | sed 's/;/\\t/' | awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$4}' | sed 's/ID=//' > genes_all.list
        Rscript ${baseDir}/bin//abstract_pop_genes.R --regionfile ${region} --genefile genes_all.list --outfile genes_abstract.list
        cat genes_abstract.list | awk 'NR>1{print \$4}' > degfile
        Rscript ${baseDir}/bin/enrich.R --degfile degfile --term2genefile ../${term2gene} --term2namefile ../${term2name} --outname ${method} --db ../${LIB}
        perl ${baseDir}/bin/extract_region.gene.eff.pl -gff ../${gff} -region ${region} -table ../${table} -out ${method}
        Rscript ${baseDir}/bin/merge_annotation.R --regionfile ${region} --genefile genes_abstract.list --gofile ../${GO_anno} --keggfile ../${KEGG_file} --outfile region_gene.txt
    """
    }
}