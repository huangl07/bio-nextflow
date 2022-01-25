#!/usr/bin/env nextflow
params.fqlist = 'fq.list'
params.outdir = "demo"
params.ref="ref.fa"
params.gff="ref.gff"
params.help = false
params.group="group.file"
params.condition="condition"
params.method="GATK"
params.anno="anno"
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --fqlist    <file>  rawdata file list
    --outdir    <dir>   output dir
    --ref   <file>  reference genome fasta file
    --gff   <file>  reference gff file 
    --group <file>  group file for different analysis
    --condition <file>  conditaion frile for different analysis
    --go_annotation <dir>  anno dir
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}
ref_file=file(params.ref)
gff_file=file(params.gff)
groups_file=file(params.group)
go_anno=file("${params.anno}/go.list");
kegg_anno=file("${params.anno}/kegg.list}");
if(params.condition){
    condition_file=file(params.condition)
}
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
            file "${name}.*.atgcn"
            file "${name}.*.qual"
        script:
        """
            fastp -i ${read1} -o ${name}_clean.1.fastq.gz -I ${read2} -O ${name}_clean.2.fastq.gz -w 8 -h ${name}.html -j ${name}.json
            perl ${baseDir}/bin/fastp.pl -i ${name}.json -o ${name}
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
        Rscript ${baseDir}/bin/ngsqc.R --base ${name}.raw.atgcn --qual ${name}.raw.qual --key ${name}.raw --od ./
        Rscript ${baseDir}/bin/ngsqc.R --base ${name}.clean.atgcn --qual ${name}.clean.qual --key ${name}.clean --od ./
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
        file "ref.chrlist" into ref_chrlist
        file "ref.fa"
        file "ref.fa" into ref_fa
    script:
    """
        ln -s ${ref} ref.fa
        samtools faidx ref.fa
        less -S ref.fa|grep ">"|grep -E "chr|chromosome|Chr|Chromosome"|cut -f 1 -d " "|sed 's/>//g' > ref.chrlist
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
process refindex {
    publishDir "${params.outdir}/03.reference" , pattern: "*"
    tag "hisat2-index"
    queue "DNA"
    cpus 16
	memory "200G"
    executor "slurm"
    input:
        val(ref) from ref_file
        val(gff) from gff_file
    output:
        file "*" into ref_index
        file "ref.gtf" into ref_gtf
    script:
    """
        gffread ${gff} -g ${ref} -T -o ref.gtf -x ref.cds.fa -y ref.pro.fa
        extract_exons.py ref.gtf > ref.exon
        extract_splice_sites.py ref.gtf > ref.ss
        hisat2-build -p 16 ${ref} --ss ref.ss --exon ref.ss ref.fa
    """
}
process snpEff {
    publishDir "${params.outdir}/03.reference/" , pattern: "ref"
    publishDir "${params.outdir}/03.reference/" , pattern: "ref.config"
    publishDir "${params.outdir}/03.reference/" , pattern: "ref.bed"
    tag "snpEff"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        val(ref) from ref_file
        val(gff) from gff_file
    output:
        file "ref.config"
        file "ref.bed"
        file "ref.bed" into ref_bed1,ref_bed2,ref_bed3
        path "ref" into snpEff_path1,snpEff_path2
        file "ref.config" into snpEff1,snpEff2
    script:
    """
        mkdir ref
        echo "data.dir= ./\nref.genome : ref" > ref.config
        ln -s ${ref} ref/sequences.fa
        ln -s ${gff} ref/genes.gff
        gffread -o ref/genes.gtf -T ${gff}
        gtf2bed ref/genes.gtf > ref.bed
        snpeff build -gtf22 -v ref -c ref.config
    """
}

process mapping {
    publishDir "${params.outdir}/04.mapping/" , pattern: "*"
    tag "mapping"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        val(ref) from ref_file
        file "*" from ref_index
        tuple val(name),file(read1),file(read2) from clean_reads
    output:
        tuple val(name),file("${name}.sam") into align_bam
        file "*.sam"
        file "*.stat" into map_stat
    script:
    """
        hisat2 -p 8 --dta -x ref.fa -1 ${read1} -2 ${read2}  -S ${name}.sam --un-conc-gz ${name}.unmapping.gz --novel-splicesite-outfile ${name}.bed --rg-id "@RG" --rg "ID:${name}" --rg "LD:${name}" --rg "PL:illumina" --rg "SM:${name}" --rg "PU:run_barcode"  --rg  "CN:MajorBio" --rg "DS:RNA" --summary-file ${name}.align.stat
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
        tuple val(name),file("${name}.sorted.bam"),file("${name}.sorted.bam.bai") into sorted_bam1,sorted_bam2,sorted_bam3,sorted_bam4,sorted_bam5,sorted_bam6,sorted_bam7,sorted_bam8,sorted_bam9
        file "${name}.sorted.bam" into sort_bams
        file "${name}.sorted.bam.bai" into sort_bams_bai
        file "${name}.*"
    script:
    """
        samtools fixmate -@ 8 -O BAM ${bam}  ${name}.fixmate.bam
        samtools sort -o ${name}.sorted.bam --output-fmt BAM -@ 8 ${bam}
        samtools index ${name}.sorted.bam
    """
}
process sampleAssemble{
    publishDir "${params.outdir}/06.sampleAssemble",pattern:"*"
    tag "sampleAssemble"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        tuple val(name),bam,bai from sorted_bam1
        file refgtf from ref_gtf
    output:
        file "${name}.primary.gtf" into assemble_gtf
    script:
    """
        stringtie -p 8 -G ${refgtf} -j 3 -c 5 -o ${name}.primary.gtf -l ${name} $bam
    """
}
process mergeAssemble{
    publishDir "${params.outdir}/07.totalmerged",pattern:"*"
    tag "totalmerged"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        file gtfs from assemble_gtf.collect()
        file refgtf from ref_gtf
    output:
        file "total.merged.gtf" into merged_gtf
    script:
        def gtf_line = gtfs.join(' ')
    """
       stringtie --merge -p 8 -G $refgtf -T 1 -c 5 -f 0.1 -o total.merged.gtf -l merged ${gtf_line}
       gffcompare -G -r $refgtf total.merged.gtf
    """
}
process reassemble{
    publishDir "${params.outdir}/08.reAssemble",pattern:"*"
    tag "reAssemble"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        tuple val(name),bam,bai from sorted_bam2
        file mergeGtf from merged_gtf
    output:
        file "${name}.final.gtf" into final_gtfs
    script:
    """
       stringtie -e -p 8 -G ${mergeGtf} -l ${name} -o ${name}.final.gtf $bam
    """
}
process prepDE{
    publishDir "${params.outdir}/09.prepDE",pattern:"*"
    tag "prepDE"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        file(finalgtfs) from final_gtfs.collect()
    output:
        file "*"
        file "gene_count_matrix.csv" into gene_count
        file "transcript_count_matrix.csv" into  trans_count
    script:
    """
        ls *.gtf|perl -ne '\$id=(split(/\\./,\$_))[0];print "\$id\\t\$_"' > gtf.list
        prepDE.py -i gtf.list --cluster --legend=legend.csv
    """
}
process mappingQC{
    publishDir "${params.outdir}/10.mappingQC",pattern:"*"
    tag "alignstat"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        file stat from map_stat.collect()
    output:
        file "*"
    script:

    """
        perl ${baseDir}/bin//mapping.pl -i ./ -o all.alignment.xls
    """
}
process geneBody_coverage{
    publishDir "${params.outdir}/10.mappingQC",pattern:"*"
    tag "geneBody"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        tuple val(name),bam,bai from sorted_bam3
        file refbed from ref_bed1
    output:
        file "*"
        file "${name}.geneBodyCoverage.txt" into coverage_file
    script:
    """
        geneBody_coverage.py -i ${bam} -r ${refbed}  -o ${name}
    """
}
process coverage_all{
    publishDir "${params.outdir}/10.mappingQC",pattern:"*"
    tag "geneBody"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        file cover_file from coverage_file.collect()
    output:
        file "*"
    script:
    """
        readlink -f *.geneBodyCoverage.txt > cover.txt
        Rscript ${baseDir}/bin/coverage.R --input cover.txt --output ${name}
    """
}
process satuation{
 publishDir "${params.outdir}/10.mappingQC",pattern:"*"
    tag "satuation"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        tuple val(name),bam,bai from sorted_bam5
        file refbed from ref_bed2
    output:
        file "*"
    script:
    """
        RPKM_saturation.py -i ${bam} -r ${refbed}  -o ${name}
        Rscript ${baseDir}/bin/saturation.R --input ${name}.eRPKM.xls --output ${name}
    """
}
process read_distribution{
 publishDir "${params.outdir}/10.mappingQC",pattern:"*"
    tag "reads_distribution"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        tuple val(name),bam,bai from sorted_bam6
        file refbed from ref_bed3
    output:
        file "*"
    script:
    """
        read_distribution.py -i ${bam} -r ${refbed} > ${name}.distribution.xls
        Rscript ${baseDir}/bin/distribution.R --input ${name}.distribution.xls --output ${name}
    """
}

process DEGene{
    publishDir "${params.outdir}/11.DEseq2",pattern:"*"
    tag "DEGene"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        file count_matrix from gene_count
        file group from groups_file
        file condition from condition_file
    output:
        path "gene"
        file "gene/*.xls" into gene_de_list
    script:
        if(condition_file){
            """
            Rscript ${baseDir}/bin/DEseq2.R --group ${group} --matrix ${count_matrix} --output gene
            """
        }else{
            """
            Rscript ${baseDir}/bin/DEseq2.R --group ${group} --matrix ${count_matrix} --output gene --contrast {condition}
            """
        }
}
process DEtrans{
    publishDir "${params.outdir}/11.DEseq2",pattern:"*"
    tag "DEtrans"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        file count_matrix from trans_count
        file group from groups_file
        file condition from condition_file
    output:
        path "transcript"
        file "transcript/*.sig.xls" into trans_de_list
    script:
        if(condition_file){
            """
            Rscript ${baseDir}/bin/DEseq2.R --group ${group} --matrix ${count_matrix} --output transcript
            """
        }else{
            """
            Rscript ${baseDir}/bin/DEseq2.R --group ${group} --matrix ${count_matrix} --output transcript --contrast {condition}
            """
        }
}
process enrichGO{
    publishDir "${params.ourdir}/19.enrichment",pattern:"*"
    tag "enrichGO-trans"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        file trans from trans_de_list
        file go from go_anno
    output:
        path "*"
    script:
        """
            Rscript ${baseDir}/bin/enrich.R -list ${trans} -annotate ${go} -output ${trans}
        """
}
process enrichKEGG{
    publishDir "${params.ourdir}/19.enrichment",pattern:"*"
    tag "enrichKEGG-trans"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        file trans from trans_de_list
        file kegg from kegg_anno
    output:
        path "*"
    script:
        """
            perl ${baseDir}/bin/kobas.pl -de ${trans} -db ${kegg} -o ${trans}.kobas
            identify.py -fg ${trans}.kobas -bg ${kegg} -o ${trans}.out
        """
}
process corPCA{
    publishDir "${params.outdir}/12.PCA-Cor",pattern:"*"
    tag "corPCA"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        file gene_matrix from gene_count
        file trans_matrix from trans_count
        file group from groups_file
    output:
        path "gene"
        path "transcript"
    script:
    """
        Rscript ${baseDir}/bin/corPCA.R --group ${group} --matrix ${gene_matrix} --output gene
        Rscript ${baseDir}/bin/corPCA.R --group ${group} --matrix ${trans_matrix} --output transcript
    """    

}
process mkdup {
    publishDir "${params.outdir}/13.mkdup/" , pattern: "*"
    tag "mkdup"
    queue "DNA"
    cpus 8
    memory "50G"
    executor "slurm"
    input:
        tuple val(name),bam,bai from sorted_bam7
        file ref from ref_fa
        file fai from ref_fai
    output:
        tuple val(name),file("${name}.mkdup.bam"),file("${name}.mkdup.bam.bai"),file("${name}.metric.txt") into mkdup_bam
        tuple val(name),file("${name}.splice.bam"),file("${name}.splice.bam.bai") into splice_bam
        file "${name}.*"
    script:
    """
        sentieon driver -t 8 -i ${bam} --algo LocusCollector --fun score_info ${name}.score.txt
        sentieon driver -t 8 -i ${bam} --algo Dedup --rmdup --score_info ${name}.score.txt --metric ${name}.metric.txt ${name}.mkdup.bam
        sentieon driver -t 8 -i ${name}.mkdup.bam -r ${ref} --algo RNASplitReadsAtJunction --reassign_mapq 255:60 ${name}.splice.bam
        #samtools markdup ${bam} ${name}.mkdup.bam -@ 8
        #samtools mkdup ${name}.mkdup.bam
        #sentieon driver -t 8 -i ${name}.mkdup.bam -r ${ref} --algo RNASplitReadsAtJunction --reassign_mapq 255:60 ${name}.splice.bam
    """
}
process bqsr{
    publishDir "${params.outdir}/14.bqsr/" , pattern: "*"
    tag "bqsr"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        tuple val(name),bam,bai from splice_bam
        file(ref) from ref_file
        file(dict) from ref_dict
        file(fai) from ref_fai
    output:
        tuple val(name),"${name}.bam","${name}.bam.bai","${name}.table" into bqsr_table
        file "*"
    script:
    """
        ln -s ${ref} ref.fa
        ln -s ${bam} ${name}.bam
        ln -s ${bai} ${name}.bam.bai
        sentieon driver -t 8 -r ref.fa  -i ${bam} --algo QualCal  ${name}.table
    """
}
process haplotyper {
    publishDir "${params.outdir}/15.haplotyper/" , pattern: "*"
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
        
        sentieon driver -t 8 -r ref.fa -i ${bam} -q ${qtable} --algo Haplotyper --ploid 2 --emit_mode GVCF ${name}.gvcf.gz        
    """
}


process gvcftyper {
    publishDir "${params.outdir}/16.gvcftyper/" , pattern: "*"
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
process snpFilterAnno{
    publishDir "${params.outdir}/17.filter/",pattern:"*"
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
    publishDir "${params.outdir}/17.filter/",pattern:"*"
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
    publishDir "${params.outdir}/17.filter/",pattern:"*"
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
process prermats{
    publishDir "${params.outdir}/18.rmats/",pattern:"*"
    tag "pre-rmats"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        file(bams) from sort_bams.collect()
        file group from groups_file
        file condition from condition_file
    output:
        file "rmats.sh" into rmats_worksh
        file "*.bam.list" into bam_list
    script:

     if(condition_file){
            """
            perl  ${baseDir}/bin/preRMATs.pl --input ./  --group ${group} --condition ${condition} -output rmats.sh -dOut ./
            """
    }else{
            """
             perl  ${baseDir}/bin/preRMATs.pl --input ./  --group ${group} -output rmats.sh -dOut ./
            """
     }

}

rmats_worksh.collectFile().splitText(by:1).set{rmats_paramter}

process runrmats{
    publishDir "${params.outdir}/18.rmats/",pattern:"*"
    tag "runramts"
    queue "DNA"
    cpus 1
    executor "slurm"
    input:
        val para from rmats_paramter
        file gtf from merged_gtf
        file bamlist from bam_list
    output:
        file "*"
    script:
    """
        run_rmats --nthread 8 --gtf ${gtf} ${para}
    """
}
