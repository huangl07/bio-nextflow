#!/usr/bin/env nextflow
params.outdir = "demo"
params.ref="ref.fa"
params.gff="ref.gff"
params.help = false
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --ref   <file>  reference genome fasta file
    --gff   <file>  reference gff file 
    --outdir   <dir>   output dir
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}
ref_file=file(params.ref)
gff_file=file(params.gff)

process fastaGenerate{
    publishDir "${params.outdir}/01.protein", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file ref from ref_file
        file gff from gff_file
    output:
        file "ref.pep.fa" into protein
        file "ref.pep.bed" into protein_bed
    script:
    """
        gffread -J -W -y ref.pep.fa ${gff} -g ${ref}
        gfftobed -m ${gff} -a transcript_id > ref.pep.bed
    """
}

protein.splitFasta(by:2000,file:true).into {fastas1;fastas2;fastas3;fastas4}

process NR {
    publishDir "${params.outdir}/02.NR", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        each fasta from fastas1
    output:
        file "*.nr.blast" into nrblast
    script:
    """
        diamond blastp --db /mnt/ilustre/users/dna/database/nr/nr.dmnd --query ${fasta} --evalue 10e-10 -k 1 --outfmt 6 qseqid sseqid salltitles evalue bitscore --threads 8 --out ${fasta.fileName}.nr.blast
    """
}
process NRanno {
    publishDir "${params.outdir}/07.Result", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file nr from nrblast.collect()
    output:
        file "NR.result" 
    script:
    """
        cat ${nr}|cut -f 1,2,3 > NR.result
    """
}
process KEGG {
    publishDir "${params.outdir}/03.KEGG", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        each fasta from fastas2
    output:
        file "*.KEGG.blast" into keggblast
    script:
    """
        diamond blastp --db /mnt/ilustre/users/dna/database/kobas/seq_pep/ko.pep.fasta.dmnd --query ${fasta} --evalue 10e-10 -k 1 --outfmt 6 --threads 8 --out ${fasta.fileName}.KEGG.blast
    """
}

process KEGGanno {
    publishDir "${params.outdir}/07.Resule", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file kegg_blast from keggblast.collect()
    output:
        file "KEGG.kobas.result" 
    script:
    """
    cat ${kegg_blast} > total.blast
    annotate.py -i total.blast -t blastout:tab -s ko -o KEGG.kobas.result -q /mnt/ilustre/users/dna/database/kobas/sqlite3/
    """
}

process Uniref {
    publishDir "${params.outdir}/04.Uniref", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        each fasta from fastas3
    output:
        file "*.uniref.blast" into unirefblast1,unirefblast2
    script:
    """
        diamond blastp --db /mnt/ilustre/users/dna/database/uniprot/uniref90.fasta.dmnd --query ${fasta} --evalue 10e-10 -k 1 --outfmt 6 qseqid sseqid salltitles evalue bitscore --threads 8 --out ${fasta.fileName}.uniref.blast
    """
}
process Unianno {
    publishDir "${params.outdir}/07.Result", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file uniref from unirefblast1.collect()
    output:
        file "Uniref.result"
    script:
    """
        cat ${uniref}|cut -f 1,2,3 > Uniref.result
    """
}
process GO {
    publishDir "${params.outdir}/05.GO", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file uniref from unirefblast2
    output:
        file "*.go" into goblast
    script:

    """
        ${baseDir}/bin/idmapping -d /mnt/ilustre/users/dna/database/idmapping_selected.tab.gz -i ${uniref} -o ${uniref}.go
    """
}
process GOanno {
    publishDir "${params.outdir}/07.Result", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file go from goblast.collect()
    output:
        file "GO.result"
    script:
    """
        cat ${go}|cut -f 1,2,3 > GO.result
    """
}
process EGGNOG { 
    publishDir "${params.outdir}/06.EGGNOG", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        each fasta from fastas4
    output:
        file "*.emapper.annotations" into eggnog
        file "*"
    script:
    """
        emapper.py -m diamond  --data_dir /mnt/ilustre/users/dna/.env/eggnog-mapper-2.1.4-2/data -o ${fasta.fileName}.eggnog -i ${fasta}
    """
}
process EGGanno { 
    publishDir "${params.outdir}/07.Result", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file eggnog from eggnog.collect()
    output:
        file "eggnog.result" 
    script:
    """
        cat ${eggnog} > eggnog.result
    """
}
