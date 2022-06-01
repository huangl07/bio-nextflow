params.outdir = "demo"

def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
   

    --gff   <file>  gff file
    --GOanno    <file>  GOanno
    --KEGGanno  <file>  KEGGanno
    --regoin    <region> regon for enrich

    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}
gff_file=file(params.gff)
GOanno=file(params.GOanno)
KEGGanno=file(params.KEGGanno)
region=file(region)
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
        file "KEGG_annotation.csv" into KEGG
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
    publishDir "${params.outdir}/04.enrich", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        tuple method,region from region
        tuple path(term2gene),path(term2name) from term_ch
        file gff from gff_file
        path LIB from EnrichDb
        file GO_anno from GOanno
        file KEGG_file from KEGG
    output:
        file "*"
    script:
    """
        mkdir ${method}
        cd ${method}
        grep mRNA ../${gff} | awk '{print \$1"\\t"\$4"\\t"\$5"\\t"\$9}' | sed 's/;/\\t/' | awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$4}' | sed 's/ID=//' > genes_all.list
        Rscript ${baseDir}/bin//abstract_pop_genes.R --regionfile ${region} --genefile genes_all.list --outfile genes_abstract.list
        cat genes_abstract.list | awk 'NR>1{print \$4}' > degfile
        Rscript ${baseDir}/bin/enrich.R --degfile degfile --term2genefile ../${term2gene} --term2namefile ../${term2name} --outname ${method} --db ../${LIB}
    """
}