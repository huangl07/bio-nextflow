#!/usr/bin/env nextflow
params.fqlist = 'fq.list'
params.outdir = "demo"
params.ref="ref.fa"
params.gff="ref.gff"
params.chrlist="chr.list"
params.help = false
params.ploid=2
params.sv = false
params.cnv = false
params.method="GATK"
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --fqlist    <file>  rawdata file list
    --outdir    <dir>   output dir
    --ref   <file>  reference genome fasta file
    --chrlist   <file>  chrlist for draw
    --gff   <file>  reference gff file 
    --ploid <file>  ploid, default 2
    --method    <str>   snp calling method default GATK,DNAscope,
    --cnv   whether do cnv default 0
    --sv    whether do sv default 0

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}
ref_file=file(params.ref)
gff_file=file(params.gff)
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
		tag "fastp"
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
    tag "fastqc"
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

process refFai {
    publishDir "${params.outdir}/03.reference" , pattern: "ref.*"
    tag "refFai"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        file(ref) from ref_file
    output:
        file "ref.fa.fai" into ref_fai
        file "ref.fa"
    script:
    """
        ln -s ${ref} ref.fa
        samtools faidx ref.fa
    """
}
process refdict {
    publishDir "${params.outdir}/03.reference" , pattern: "*"
    tag "refdict"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        val(ref) from ref_file
    output:
        file "ref.dict" into ref_dict
        file "ref.stat" into ref_stat
        file "ref.stat"
        file "ref.dict"
    script:
    """
        ln -s ${ref} ref.fa
        samtools dict ${ref} -o ref.dict
        perl ${baseDir}/bin/refstat.pl -ref ref.fa -output ref.stat
    """
}
process bwaindex {
    publishDir "${params.outdir}/03.reference" , pattern: "*.{0123,ann,amb,64,pac}"
    tag "bwa-index"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        val(ref) from ref_file
    output:
        file "*.{0123,ann,amb,64,pac}" into bwa_index
        file "*.{0123,ann,amb,64,pac}" 
    script:
    """
        ln -s ${ref} ref.fa
        bwa-mem2.avx index ref.fa
    """
}
process snpEff {
    publishDir "${params.outdir}/03.reference/" , pattern: "ref"
    publishDir "${params.outdir}/03.reference/" , pattern: "ref.config"
    tag "snpEff"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        val(ref) from ref_file
        val(gff) from gff_file
    output:
        file "ref.config"
        path "ref" into snpEff_path1,snpEff_path2
        file "ref.config" into snpEff1,snpEff2
    script:
    """
        mkdir ref
        echo "data.dir= ./\nref.genome : ref" > ref.config
        ln -s ${ref} ref/sequences.fa
        ln -s ${gff} ref/genes.gff
        gffread -o ref/genes.gtf -T ${gff}
        snpeff build -gtf22 -v ref -c ref.config
    """
}

process mapping {
    publishDir "${params.outdir}/04.mapping/" , pattern: "*.bam"
    tag "mapping"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        val(ref) from ref_file
        file "*" from bwa_index
        tuple val(name),file(read1),file(read2) from clean_reads
    output:
        tuple val(name),file("${name}.bam") into align_bam
        file "*.bam"
    script:
    """
        ln -s ${ref} ref.fa
        bwa-mem2.avx mem -M -a -t 8 -R "@RG\\tID:${name}\\tLG:${name}\\tLB:2\\tPL:illumina\\tSM:${name}\\tPU:run_barcode\\tCN:MajorBio\\tDS:reseq" ref.fa ${read1} ${read2}| samtools view -bS - > ${name}.bam
    """
}

process bamsort {
    publishDir "${params.outdir}/05.sort/" , pattern: "*"
    tag "bamsort"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        tuple val(name),bam from align_bam
    output:
        tuple val(name),file("${name}.sorted.bam"),file("${name}.sorted.bam.bai") into sorted_bam
        file "${name}.*"
    script:
    """
        samtools sort -o ${name}.sorted.bam --output-fmt BAM -@ 8 ${bam}
        samtools index ${name}.sorted.bam
    """
}
process mkdup {
    publishDir "${params.outdir}/06.mkdup/" , pattern: "*"
    tag "mkdup"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        tuple val(name),bam,bai from sorted_bam
    output:
        tuple val(name),file("${name}.mkdup.bam"),file("${name}.mkdup.bam.bai"),file("${name}.metric.txt") into mkdup_bam
        file "${name}.*"
    script:
    """
        sentieon driver -t 8 -i ${bam} --algo LocusCollector --fun score_info ${name}.score.txt
        sentieon driver -t 8 -i ${bam} --algo Dedup --rmdup --score_info ${name}.score.txt --metric ${name}.metric.txt ${name}.mkdup.bam
    """
}
process realign {
    publishDir "${params.outdir}/07.realign/" , pattern: "*"
    tag "realign"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        file "ref.fa.fai" from ref_fai
        file "*" from bwa_index
        file(ref) from ref_file
        tuple val(name),bam,bai,metric from mkdup_bam
    output:
        tuple val(name),file("${name}.realign.cram"),file("${name}.realign.cram.crai"),file("${name}.metric") into realign_cram1,realign_cram3,realign_cram4,realign_cram5
        tuple val(name),file("${name}.realign.bam"),file("${name}.realign.bam.bai"),file("${name}.metric") into realign_cram2
        file "${name}.realign.cram" into realign_bam
        file "${name}.realign.cram.crai" into realign_index
        file "${name}.*"
    script:
    """
        ln -s ${ref} ref.fa
        ln -s ${metric} ${name}.metric
        sentieon driver -t 8 -i ${bam} -r ref.fa --temp_dir ./ --algo Realigner ${name}.realign.bam
        samtools view -T ref.fa -C -o ${name}.realign.cram ${name}.realign.bam
        samtools index ${name}.realign.cram
    """
}
process insertcoverage {
    publishDir "${params.outdir}/08.mapstat/" , pattern: "*"
    tag "insert_coverage"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        tuple val(name),bam,bai,mertric from realign_cram1
        val(dict) from ref_dict
    output:
        file "*"
    script:
    """
        samtools stats -@ 8 $bam > ${name}.all.mapstat 
        perl ${baseDir}/bin//map_stat.pl -b ${name}.all.mapstat -m ${mertric} -i ${name}.insert -c ${name}.coverage -d ${dict} -o ${name}.result.stat -k ${name}
        Rscript ${baseDir}/bin/insert_size.R --i ${name}.insert --o ${name}.insert
        Rscript ${baseDir}/bin/coverage_depth.R --i ${name}.coverage --o ${name}.depth
    """
}
process depth_stat {
    publishDir "${params.outdir}/08.mapstat/" , pattern: "*"
    tag "depth"
    queue "DNA"
    cpus 8
    executor "slurm"
    memory "100G"
    input:
        tuple val(name),bam,bai,mertric from realign_cram2
        val(dict) from ref_dict
        val(chrlist) from chr_list
    output:
        file "*"
    script:
    """
        sambamba depth window -w 100000 -t 8 -o ${name}.stat $bam 
        Rscript ${baseDir}/bin/genomeCoveragehorizontalArea.R --infile ${name}.stat  --idfile ${chrlist} --outfile ${name}.genome.coverage --group.col 1 --x.col 2 --y.col 4 --x.lab Sequence-Position --y.lab AverageDepth-log2 --skip 0 --unit kb --log2
    """
}
if(params.method == "GATK"){
    process bqsr{
        publishDir "${params.outdir}/09.bqsr/" , pattern: "*"
        tag "bqsr"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            tuple val(name),bam,bai,mertric from realign_cram3
            file(ref) from ref_file
            file(dict) from ref_dict
            file(fai) from ref_fai
        output:
            tuple val(name),"${name}.cram","${name}.cram.crai","${name}.table" into bqsr_table
            file "*"
        script:
        """
            ln -s ${ref} ref.fa
            ln -s ${bam} ${name}.cram
            ln -s ${bai} ${name}.cram.crai
            sentieon driver -t 8 -r ref.fa  -i ${bam} --algo QualCal  ${name}.table
        """
    }
    process haplotyper {
        publishDir "${params.outdir}/09.haplotyper/" , pattern: "*"
        tag "haplotyper"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            tuple val(name),bam,bai,qtable from bqsr_table
            val(ref) from ref_file
            file(dict) from ref_dict
            file(fai) from ref_fai
        output:
            file "${name}.gvcf.gz" into gvcf
            file "${name}.gvcf.gz.tbi" into gvcf_tbi
            file "*"
        script:
        """
            ln -s ${ref} ref.fa
            
            sentieon driver -t 8 -r ref.fa -i ${bam} -q ${qtable} --algo Haplotyper --ploid ${params.ploid} --emit_mode GVCF ${name}.gvcf.gz        
        """
    }   
    process gvcftyper {
        publishDir "${params.outdir}/10.gvcftyper/" , pattern: "*"
        tag "gvcftype"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            file(gvcfs) from gvcf.collect()
            file(tbis) from gvcf_tbi.collect()
            file(ref) from ref_file
            file(dict) from ref_dict
            file(fai) from ref_fai
        output:
            file "pop.variant.vcf" into vcf
            file "pop.variant.vcf.idx" into vcf_index
            file "*"
        script:
            def gvcf_line = gvcfs.collect{ "-v $it" }.join(' ')
        """
            ln -s ${ref} ref.fa
            sentieon driver -t 8 -r ref.fa --algo GVCFtyper ${gvcf_line}  pop.variant.vcf
        """
    }
}else{
    process DNAscope {
        publishDir "${params.outdir}/09.DNAscope/" , pattern: "*"
        tag "DNAscope"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            tuple val(name),bam,bai,qtable from bqsr_table
            val(ref) from ref_file
            file(dict) from ref_dict
            file(fai) from ref_fai
        output:
            file "${name}.vcf.gz" into vcf_list
            file "${name}.vcf.gz.tbi" into gvcf_tbi
            file "*"
        script:
        """
            ln -s ${ref} ref.fa
            sentieon driver -t 8 -r ref.fa -i ${bam} --algo DNAscope --ploid ${params.ploid}  ${name}.vcf.gz       
        """
    }
    process Vcfmerge {
        publishDir "${params.outdir}/10.Vcfmerge/" , pattern: "*"
        tag "Vcfmerge"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            file(gvcfs) from gvcf.collect()
            file(tbis) from gvcf_tbi.collect()
            file(ref) from ref_file
            file(dict) from ref_dict
            file(fai) from ref_fai
        output:
            file "pop.variant.vcf.gz" into vcf
            file "pop.variant.vcf.gz.idx" into vcf_index
            file "*"
        script:
            def vcfs = vcf_list.collect{ "$it" }.join(' ')
        """
            ln -s ${ref} ref.fa
            bcftools merge -o pop.variant.vcf.gz -O z  --threads 8 --missing-to-ref $vcfs
            bcftools index pop.variant.vcf.gz
        """
    }
}

 
if(params.sv){
    process bamPEforManta {
        publishDir "${params.outdir}/14.bamPE/" , pattern: "*"
        tag "bamPEforManta"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            tuple val(name),bam,bai,mertric from realign_cram4
        output:
            file "${name}.PE.cram"  into PEbam
            file "${name}.PE.cram.crai" into PEindex
            file "${name}.PE.cram"
            file "${name}.PE.cram.crai"
        script:
        """
            samtools view -b --output-fmt CRAM -f 2 ${bam} >${name}.PE.cram
            samtools index ${name}.PE.cram
        """
    }
    process runManta {
        publishDir "${params.outdir}/15.svCalling/" , pattern: "*"
        tag "runManta"
        queue "DNA"
        cpus 16
        executor "slurm"
        input:
            file(bam) from PEbam.collect()
            file(bai) from PEindex.collect()
            file(ref) from ref_file
            file(dict) from ref_dict
            file(fai) from ref_fai
        output:
            file "Manta/results/variants/diploidSV.vcf.gz" into sv_vcf
            file "Manta/results/variants/diploidSV.vcf.gz.tbi" into sv_index
            file "sv.xls"
            file "Manta/results/variants/diploidSV.vcf.gz"
            file "Manta/results/variants/diploidSV.vcf.gz.tbi"
        script:
            def bams = bam.collect{"--bam=$it" }.join(' ')
        """
            ln -s ${ref} ref.fa
            configManta.py --referenceFasta=ref.fa ${bams} --runDir Manta
            Manta/runWorkflow.py -j 16 -m local 
            perl ${baseDir}/bin/sv-stat.pl -i Manta/results/variants/diploidSV.vcf.gz -o sv.xls
        """
    }
}
if(params.cnv){
    process CreatePon {
        publishDir "${params.outdir}/16.createpon/" , pattern: "*"
        tag "CreatePon"
        queue "DNA"
        cpus 16
        executor "slurm"
        input:
            file(bams) from realign_bam.collect()
            file(indexs) from realign_index.collect()
            file(ref) from ref_file
            file(dict) from ref_dict
            file(fai) from ref_fai
        output:
            file "total.pon" into pon_file
            file "*"
        script:
            def bams = bams.collect{"-i $it" }.join(' ')
        """
            ln -s ${ref} ref.fa
            less -S ref.fa.fai|perl -ne 'chomp;@a=split;print "\$a[0]\\t1\\t\$a[1]\\t\$a[0]\\n";' > target_bed
            sentieon driver -t 16 -r ref.fa  ${bams} --algo CNV --target target_bed --create_coverage total.coverage
            sentieon driver -t 16 -r ref.fa  ${bams} --algo CNV --coverage total.coverage --create_pon total.pon
        """
    }
    process CNVcalling {
        publishDir "${params.outdir}/17.cnvcalling/" , pattern: "*"
        tag "CNVcalling"
        queue "DNA"
        cpus 16
        executor "slurm"
        input:
            tuple val(name),bam,bai,mertric from realign_cram5
            file(ref) from ref_file
            file(dict) from ref_dict
            file(fai) from ref_fai
            file(pon) from pon_file
        output:
            file "*"
            file "*.cnv" into cnvfiles
        script:
        """
        ln -s ${ref} ref.fa
        sentieon driver -t 16 -r ref.fa -i ${bam} --algo CNV --pon ${pon} ${name}.cnv
        perl  ${baseDir}/bin/cnv-stat.pl -input ${name}.cnv -output ${name}.cnv.xls
        """
    }
}
process snpFilterAnno{
    publishDir "${params.outdir}/18.filter/",pattern:"*"
    tag "snpfilter"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        file vcf from  vcf
        file vcfindex from vcf_index
        path ref_path from snpEff_path1
        file ref_config from snpEff1
    output:
        file "pop.snp.filtered.vcf.gz" into filter_snp_vcf
        file "pop.snp.filtered.anno.vcf" into anno_snp_vcf
        file "snpEff.snp.csv" 
        file "snp.sample.xls";
        file "snp.effects.xls";
        file "snp.region.xls";
        file "snp.summary.html";
    script:
    """
        gatk SelectVariants -V ${vcf}  -select-type SNP  -O pop.snp.vcf.gz --select-type-to-exclude SYMBOLIC 
        gatk VariantFiltration -V pop.snp.vcf.gz  --filter-expression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0"  --filter-name "HardFiltered" --set-filtered-genotype-to-no-call True -O pop.snp.filtered.vcf.gz
        snpeff eff -v ref -htmlStats snp.summary.html -csvStats snpEff.snp.csv  -c ${ref_config} pop.snp.filtered.vcf.gz > pop.snp.filtered.anno.vcf
        perl  ${baseDir}/bin/snpEff.pl -i snpEff.snp.csv -output snp
    """
}
process indelFilterAnno{
    publishDir "${params.outdir}/18.filter/",pattern:"*"
    tag "indelfilter"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        file vcf from  vcf
        file vcfindex from vcf_index
        path ref_path from snpEff_path2
        file ref_config from snpEff2
    output:
        file "pop.indel.filtered.vcf.gz" into filter_indel_vcf
        file "pop.indel.filtered.anno.vcf" into anno_indel_vcf
        file "snpEff.indel.csv" 
        file "indel.sample.xls" 
        file "indel.effects.xls";
        file "indel.region.xls";
        file "indel.summary.html";
    script:
    """
        gatk SelectVariants -V ${vcf}  -select-type INDEL  -O pop.indel.vcf.gz --select-type-to-exclude SYMBOLIC 
        gatk VariantFiltration -V pop.indel.vcf.gz  --filter-expression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0"  --filter-name "HardFiltered" --set-filtered-genotype-to-no-call True -O pop.indel.filtered.vcf.gz
        snpeff eff -v ref -htmlStats indel.summary.html -csvStats snpEff.indel.csv  -c ${ref_config} pop.indel.filtered.vcf.gz > pop.indel.filtered.anno.vcf
        perl  ${baseDir}/bin/snpEff.pl -i snpEff.indel.csv -output indel
    """    
}
process mergeVcfs{
    publishDir "${params.outdir}/18.filter/",pattern:"*"
    tag "indelfilter"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        file indel from  anno_indel_vcf
        file snp   from  anno_snp_vcf
    output:
        file "pop.final.vcf.gz" into final_vcf
    script:
    """
        #gatk IndexFeatureFile -I ${indel}
        #gatk IndexFeatureFile -I ${snp}
        gatk MergeVcfs  -I ${indel} -I ${snp} -O pop.final.vcf.gz 
    """    
}



