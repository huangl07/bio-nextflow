#!/usr/bin/env nextflow
params.outdir = "demo"
params.ref="ref.fa"
params.gff="ref.gff"
params.help = false
params.method="GATK"
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --outdir    <dir>   output dir
    --query     <file> query fasta
    --dblist    <file> dblist in a file for blast
    --remake    <boolean> default false
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

Channel.from(file(params.query)).splitFasta(by:2000,file:true).set{fastas}
Channel.from(file(params.dblist)).splitText(by:1).set{dbs}



process NR{
    publishDir "${params.outdir}/02.nrblast", pattern:"*"
    tag "NR"
    queue "DNA"
    cpus 16
    executor "slurm"
    memory "100G"
    input:
        each fasta from fastas
        each db from dbs
    output:
        file "*"
        file "nr.blast" into nr_result
    script:
        
    """
        blastp  -num_threads 16 -max_target_seqs 1 -out nr.blast -evalue 1e-5 -query ${fasta} -outfmt 6 -db ${db}
    """
}

nr_result.collectFile(name:"blast.out")