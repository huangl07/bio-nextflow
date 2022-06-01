#!/usr/bin/env nextflow
params.outdir = "demo"
params.help = false
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --mark    <file>  input vcf file
    --map     <file>   output dir
    --trt     <file>  input group file    
    --outdir  <dir>   output dir
    --popt  <str>   population type


    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}
 mark=file(params.mark)
 trt =file(params.trt)
 map =file(params.map)
 //gff_file=file(params.gff)
 //GOanno=file(params.GOanno)
 //KEGGanno=file(params.KEGGanno)

//process Orgdb {
//    publishDir "${params.outdir}/00.orgdb", pattern:"*"
//    executor 'slurm'
//    queue 'DNA'
//    cpus 2
//    cache 'lenient'
//    input:
//        file GO_anno from GOanno
//        file KEGG_anno from KEGGanno
//    output:
//        tuple path("TERM2GENE.txt"),path("TERM2NAME.txt") into term_ch
//        path(LIB) into EnrichDb
//        file "KEGG_annotation.csv" into KEGG
//    script:
//        """
//        mkdir LIB
//        awk -v RS='----+' 'NR==2' ${KEGGanno} | awk -v RS="////" '{print}' | awk NF > KEGG_raw.result
//        python  ${baseDir}/bin//KEGG_rearrange.py -filename KEGG_raw.result -outfile KEGG_annotation.csv
//        cat KEGG_annotation.csv | awk -F "\\t" 'NR>1{print \$4"\\t"\$1}' | grep -v "^-" > TERM2GENE.txt
//        cat KEGG_annotation.csv | awk -F "\\t" 'NR>1{print \$4"\\t"\$3}' | grep -v "^-" > TERM2NAME.txt
//        Rscript ${baseDir}/bin//make_orgdb.R --gofile ${GOanno} --keggfile KEGG_annotation.csv --outpath .
//        R CMD INSTALL org.Ddemo.eg.db --library=./LIB/
//        """
//} 
process dataprepair {
    publishDir "${params.outdir}/01.dataprepair", pattern:"*"
    executor 'slurm'
    queue 'DNA'
    cpus 2
    cache 'lenient'
    input:
        file mark from mark
        file map from map
        file trt from trt

    output:
        file "total.rqtl.csv" into mark_file1,mark_file2
        file "total.map" into  map_file1,map_file2
        file "total.qtl.list" into qtl
        file "total.btl.list" into btl
    script:
    if(params.popt == "CP"){
        """
        ln -s ${map} total.map
        perl ${baseDir}/bin/map2rqtlCP.pl -l $mark -o total.rqtl.csv
        Rscript ${baseDir}/bin/trt.R -i $trt -o tmp
        perl ${baseDir}/bin/trt.pl -i tmp.qtl.txt -o total.qtl --popt ${params.popt}
        perl ${baseDir}/bin/trt.pl -i tmp.btl.txt -o total.btl --popt ${params.popt}
        """
    }else{
        """
        perl ${baseDir}/bin/map2rqtl.pl -l $mark -m $map -o total.rqtl.csv
        ln -s ${map} total.map
        Rscript ${baseDir}/bin/trt.R -i $trt -o tmp
        perl ${baseDir}/bin/trt.pl -i tmp.qtl.txt -o total.qtl --popt ${params.popt}
        perl ${baseDir}/bin/trt.pl -i tmp.btl.txt -o total.btl --popt ${params.popt}
        """
    }
} 

qtl.splitCsv(header:false,sep:'\t').groupTuple().set{qtl_list}
btl.splitCsv(header:false,sep:'\t').groupTuple().set{btl_list}
process qtl {
    publishDir "${params.outdir}/02.qtl", pattern:"*"
    executor 'slurm'
    queue 'DNA'
    cpus 2
    cache 'lenient'
    input:
        file mark from mark_file1
        file map from  map_file1
        tuple qname,qtl_file from qtl_list
    output:
        path "${qname}"
        tuple qname,"${qname}/${qname}.qtl-result.result"  into qtl_result
    script:
    if(params.popt == "CP"){
        """
        Rscript ${baseDir}/bin/qtl-CP.R --map ${map} --loc ${mark} --trt ${qtl_file[0]} --out ${qname}
        """
    }else{
        """
        Rscript ${baseDir}/bin/qtl-NOCP.R --map ${map} --loc ${mark} --trt ${qtl_file[0]} --out ${qname}
        """
    }
} 
process btl {
    publishDir "${params.outdir}/03.btl", pattern:"*"
    executor 'slurm'
    queue 'DNA'
    cpus 2
    cache 'lenient'
    input:
        file "total.rqtl.csv" from mark_file2
        file "total.map" from  map_file2
        tuple bname,btl_file from btl_list
    output:
        path "${qname}"
        tuple qname,"${qname}/${qname}.btl-result.result"  into btl_result
    script:
   if(params.popt == "CP"){
        """
        Rscript ${baseDir}/bin/qtl-CP.R --map ${map} --loc ${mark} --trt ${btl_file[0]} --out ${qname} --btl 1 
        """
    }else{
        """
        Rscript ${baseDir}/bin/qtl-NOCP.R --map ${map} --loc ${mark} --trt ${btl_file[0]} --out ${qname}  --btl 1 
        """
    }
} 

//hannel.from(qtl_result).mix(btl_result).set{regions}
//rocess enrich{
//   publishDir "${params.outdir}/01.dataprepair", pattern:"*"
//   executor 'slurm'
//   queue 'DNA'
//   cpus 2
//   input:
//       tuple qname,file(result) from regions
//       path LIB from EnrichDb
//       file GO_anno from GOanno
//       file table from table_file
//       file KEGG_file from KEGG
//       file gff from gff_file
//   output:
//       file "*"
//   script:
//   """
//   mkdir ${method}
//   cd ${method}
//   grep -v lod ${result} | perl -ne '@a=split;@b=split("_",\$a[5]);@c=split("_",\$a[6]);print join("\\t",\$b[0],sort(\$b[1],\$c[1])),"\\n"' > ${qname}.region
//   grep mRNA ../${gff} | awk '{print \$1"\\t"\$4"\\t"\$5"\\t"\$9}' | sed 's/;/\\t/' | awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$4}' | sed 's/ID=//' > genes_all.list
//   Rscript ${baseDir}/bin//abstract_pop_genes.R --regionfile ${qname}.region --genefile genes_all.list --outfile genes_abstract.list
//   cat genes_abstract.list | awk 'NR>1{print \$4}' > degfile
//   Rscript ${baseDir}/bin/enrich.R --degfile degfile --term2genefile ../${term2gene} --term2namefile ../${term2name} --outname ${method} --db ../${LIB}
//   perl ${baseDir}/bin/extract_region.gene.eff.pl -gff ../${gff} -region ${qname}.region -table ../${table} -out ${method}
//   Rscript ${params.scripts}/merge_annotation.R --regionfile ${qname}.region --genefile genes_abstract.list --gofile ../${GO_anno} --keggfile ../${KEGG_file} --outfile region_gene.txt
//   """
//



workflow.onComplete {
    println "Pipeline completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
