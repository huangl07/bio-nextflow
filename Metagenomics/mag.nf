#!/usr/bin/env nextflow
params.fqlist = 'fq.list'
params.outdir = "demo"
params.ref="None"
params.chrlist="chr.list"
params.help = false
params.ploid=2
params.sv = false
params.cnv = false
params.method="megahit"
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --fqlist    <file>  rawdata file list
    --outdir    <dir>   output dir
    --ref       <file>  reference host genome fasta 
    --method    <str>   megahit or spades


    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}
ref_file=file(params.ref)
chr_list=file(params.chrlist)
reads_count = Channel.from(file(params.fqlist))
              .splitCsv(header:false,sep:'\t')
              .groupTuple()


process mergefq{
    publishDir "${params.outdir}/01.mergeFQ", pattern:"*.gz"
    queue "DNA"
    executor "slurm"
    input:
        tuple val(name),read1,read2 from reads_count
    output:
        tuple val(name),"${name}.R1.fastq.gz","${name}.R2.fastq.gz" into reads
        file "${name}*.fastq.gz"
    script:
        fqsize=read1.size() 
        fq1=read1.join(" ") 
        fq2=read2.join(" ") 
    if (fqsize > 1)
        """
            cat ${fq1} > ${name}.R1.fastq.gz
            cat ${fq2} > ${name}.R2.fastq.gz
        """
    else
        """
            ln -s ${read1[0]} ${name}.R1.fastq.gz
            ln -s ${read2[0]} ${name}.R2.fastq.gz
        """
}

//return  1
// 

process fastp {
	    publishDir "${params.outdir}/02.cleanFQ" , pattern: "*clean*.gz"
        publishDir "${params.outdir}/02.cleanFQ/QC" , pattern: "*.html"
        publishDir "${params.outdir}/02.cleanFQ/QC" , pattern: "*.json"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            tuple val(name), read1,read2 from reads
        output:
            tuple val(name),"${name}_clean.1.fastq.gz","${name}_clean.2.fastq.gz" into clean_reads
            tuple val(name),file("${name}.json")  into json
            file "${name}.json"
            file "${name}.html"
        script:
        """
            fastp -i ${read1} -o ${name}_clean.1.fastq.gz -I ${read2} -O ${name}_clean.2.fastq.gz -w 8 -h ${name}.html -j ${name}.json
        """
}


process fastqc {
    publishDir "${params.outdir}/02.cleanFQ/QC" , pattern: "*"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        tuple val(name), file(jsonfile) from json
    output:
        file "${name}.*.atgcn"
        file "${name}.*.qual"
        file "${name}.*.pdf"
        file "${name}.*.png"
        file "${name}.stat"
    script:
    """
        perl ${baseDir}/bin/fastp.pl -i ${jsonfile} -o ${name}
        Rscript ${baseDir}/bin/ngsqc.r --base ${name}.raw.atgcn --qual ${name}.raw.qual --key ${name}.raw --od ./
        Rscript ${baseDir}/bin/ngsqc.r --base ${name}.clean.atgcn --qual ${name}.clean.qual --key ${name}.clean --od ./
    """
}

if(params.ref){
    process refFai {
        publishDir "${params.outdir}/03.reference" , pattern: "ref.*"
        queue "DNA"
        cpus 1
        executor "slurm"
        input:
            file(ref) from ref_file
        output:
            file "ref.*" into ref_fai
            file "ref.fa"
        script:
        """
            ln -s ${ref} ref.fa
            bowtie2-build --threads 8 ref.fa ref
        """
    }
    process mapping {
        publishDir "${params.outdir}/04.mapping/" , pattern: "*.bam"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            val(ref) from ref_file
            file "*" from ref_fai
            tuple val(name),file(read1),file(read2) from clean_reads
        output:
            tuple val(name),file("${name}.unmap.1.fastq"),file("${name}.unmap.2.fastq"),file("${name}.unmap.single.fastq") into assemble_reads1,assemble_reads2,assemble_reads3,assemble_reads4,assemble_reads5
            file "*"
        script:
        """
            bowtie2 -x ref -1 ${read1} -2 ${read2} -S ${name}.sam  2>${name}.mapping.log 
            samtools fastq  -@ 8 -f 4 ${name}.sam -1 ${name}.unmap.1.fastq -2 ${name}.unmap.2.fastq -s ${name}.unmap.single.fastq
        """
    }
}else{
    clean_reads.combine("nonfastq").into{assemble_reads1;assemble_reads2;assemble_reads3;assemble_reads4}
}

if(params.method == "megahit"){
    process assemble_megahit {
        publishDir "${params.outdir}/06.megahit/" , pattern: "*"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            tuple val(name),file(read1),file(read2),reads3 from assemble_reads1
        output:
            tuple val(name),"megahit_out_${name}/${name}.contigs.fa" into assemble_fasta
            file "*"
            file "megahit_out_${name}/${name}.contigs.fa"
        script:
        if(params.ref == "None"){
            read3=file(${reads3})
            """
            megahit -1 ${read1} -2 ${read2} -r ${read3} -o megahit_out_${name} --out-prefix ${name} -t 8
            sed -i \'s/>/>${name}_/g\' megahit_out_${name}/${name}.contigs.fa
            """
        }else{
            """
            megahit -1 ${read1} -2 ${read2} -o megahit_out_${name} --out-prefix ${name} -t 8
            sed -i \'s/>/>${name}_/g\' megahit_out_${name}/${name}.contigs.fa
            """
        }
    }
}else{
    process assemble_spades {
        publishDir "${params.outdir}/07.spades/" , pattern: "*"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            tuple val(name),file(read1),file(read2),reads3 from assemble_reads2
        output:
            tuple val(name),file("spades_out_${name}/${name}.scaffolds.fasta")  into assemble_fasta
            file "*"
            file "spades_out_${name}/${name}.scaffolds.fasta"
        script:
        if(params.ref == "None"){
            read3=file(reads3)
            """
            spades.py --meta -1 ${read1} -2 ${read2} -r ${read3} -o spades_out_${name} -t 8
            cp spades_out_${name}/scaffolds.fasta spades_out_${name}/${name}.scaffolds.fasta
            sed -i \'s/>/>${name}_/g\' spades_out_${name}/scaffolds.fasta
            """
        }else{

            """
            spades.py --meta  -1 ${read1} -2 ${read2} -o spades_out_${name}  -t 8
            cp spades_out_${name}/scaffolds.fasta spades_out_${name}/${name}.scaffolds.fasta
            sed -i \'s/>/>${name}_/g\' spades_out_${name}/scaffolds.fasta
            """
        }
    }
}
krakendb="/mnt/ilustre/users/dna/database/kraken/kraken-database/k2_standard/"
process kranken {
    publishDir "${params.outdir}/08.kranken/" , pattern: "*"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        tuple val(name),file(read1),file(read2),reads3 from assemble_reads3
    output:
        file "*" into kranken
        file "${name}.*"
    script:

    if(params.ref == "None"){
        read3=file(reads3)
        """
        cat ${read1} ${read2} ${read3}> ${name}.fastq
        kraken2 -db ${krakendb}  --threads 8  ${name}.fastq --report ${name}.report > ${name}.out
        bracken -d ${krakendb}  -i  ${name}.report -o  ${name}.bracken -t 8
        """

    }else{
        """
        cat ${read1} ${read2} > ${name}.fastq
        kraken2 -db ${krakendb}  --threads 8 ${name}.fastq --report ${name}.report > ${name}.out
        bracken -d ${krakendb}  -i ${name}.report -o ${name}.bracken -t 8
        """
    }
}

mpadb="/mnt/ilustre/users/dna/database/mpa-database/"
process metaphalan {
    publishDir "${params.outdir}/09.metaphalan/" , pattern: "*"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        tuple val(name),file(read1),file(read2),reads3 from assemble_reads4
    output:
        file "*" into metaphalan
        file "${name}.*"
    script:

    if(params.ref == "None"){
        read3=file(reads3)
        """
        cat ${read1} ${read2} ${read3}> ${name}.fastq
        metaphlan --input_type fastq ${name}.fastq ${name}.metaphalan --bowtie2db ${mpadb} --nproc 8
        """
    }else{
        """
        cat ${read1} ${read2} > ${name}.fastq
        metaphlan --input_type fastq  ${name}.fastq ${name}.metaphalan --bowtie2db ${mpadb} --nproc 8
        """
    }
}
process metaquast {
    publishDir "${params.outdir}/10.binning/" , pattern: "*"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        tuple val(name),file(fasta)  from assemble_fasta.collect()
    output:
        file "*" 
    script:
        def fastas=fasta..collect().join(" ") 
    """
    metaquast.py -m 2000 -t 8  ${fastas} -o metaquast
    """
}
process cluster {
    publishDir "${params.outdir}/10.binning/" , pattern: "*"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        tuple val(name),file(fasta) from assemble_fasta.collect()
    output:
        file "*" into cluster_fasta
    script:
    def fastas=fasta..collect().join(" ")
    """
    mmseqs easy-cluster *.fa --threads 8 mmseqs
    """
}


process binning {
    publishDir "${params.outdir}/10.binning/" , pattern: "*"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
    output:
        file "*" into binning
    script:
    """
    touch 1.txt
    """
}