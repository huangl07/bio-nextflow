#!/usr/bin/env nextflow
params.outdir = "demo"
params.help = false
params.method="mutect2"
params.exome=false
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --normal   <file>  normal bam file
    --tumor    <file>  tumor  bam file
    --out   <dir>   output dir
    --ref   <path>  reference genome path
    --method    <str>   mutect2 or dnascope
    --exome     <bool>  exome or wgs
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}
dbsnp="/mnt/ilustre/users/dna/database/dbsnp/GCF_000001405.39.gz"
cosmic="/mnt/ilustre/users/dna/snpdatabase/CosmicCodingMuts.vcf.gz"
ref_ch1=channel.fromPath(params.ref)
ref_ch2=channel.fromPath(params.ref)
ref_ch3=channel.fromPath(params.ref)
ref_ch4=channel.fromPath(params.ref)
ref_ch5=channel.fromPath(params.ref)
normal=file(params.normal)
tumor=file(params.tumor)
channel.from([["normal",normal],["tumor",tumor]]).combine(ref_ch1).into{bamlist1;bamlist2;bamlist3;bamlist4}

process bqsr{
     publishDir "${params.outdir}/01.bqsr/" , pattern: "*"
     tag "bqsr"
     queue "DNA"
     cpus 8
     executor "slurm"
     input:
        tuple type,file(bam),ref from bamlist1
        val dbsnp from dbsnp
     output:
        file "*" into bqsr
     script:

     """ 
        samtools index ${bam}
        sentieon driver -t 8 -r ${ref}/ref.fa  -i ${bam} --algo QualCal  ${bam.simpleName}.table --known_sites $dbsnp
     """
}
if(params.method =="mutect2"){
    process mutect2 {
        publishDir "${params.outdir}/02.Mutect2/" , pattern: "*"
        tag "gvcftype"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            file bqsr from bqsr.collect()
            file normal from normal
            file tumor from tumor
            path ref from ref_ch2
        output:
            file "somatic.filter.vcf.gz" into vcf
            file "somatic.filter.vcf.gz.tbi" into vcfindex
            file "*"
        script:
        """
            sentieon driver -t 8 -i ${normal} -q ${normal.simpleName}.table -i ${tumor}  -q ${tumor.simpleName}.table -r ${ref}/ref.fa --algo TNhaplotyper2  somatic.vcf --normal_sample ${normal.simpleName} --tumor_sample ${tumor.simpleName}
            sentieon driver -r 03.reference/ref.fa --algo TNfilter --normal_sample ${normal.simpleName} --tumor_sample ${tumor.simpleName} -v somatic.vcf somatic.filter.vcf
            bcftools view -f PASS -O z somatic.filter.vcf  > somatic.filter.vcf.gz
            tabix somatic.filter.vcf.gz
        """
    }
}else{
        process TNScope {
        publishDir "${params.outdir}/02.Mutect2/" , pattern: "*"
        tag "gvcftype"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            file bqsr from bqsr.collect()
            file normal from normal
            file tumor from tumor
            path ref from ref_ch2
        output:
            file "somatic.filter.vcf.gz" into vcf
            file "somatic.filter.vcf.gz.tbi" into vcfindex
            file "*"
        script:
        """
            sentieon driver -t 8 -i ${normal} -q ${normal.simpleName}.table -i ${tumor}  -q ${tumor.simpleName}.table -r ${ref}/ref.fa --algo TNScope  somatic.vcf --normal_sample ${normal.simpleName} --tumor_sample ${tumor.simpleName} --dbsnp ${dbsnp} --cosmic ${cosmic}
            bcftools view -f PASS -O z somatic.vcf  > somatic.filter.vcf.gz
            tabix somatic.filter.vcf.gz
        """
    }
}

process snpeff{
    publishDir "${params.outdir}/03.snpeff/",pattern:"*"
    tag "snpfilter"
    queue "DNA"
    cpus 8
    executor "slurm"
    input:
        file vcf from  vcf
        file index from vcfindex
        file ref from ref_ch3
    output:
        file "*"
    script:
    """
        snpeff eff -v ref -htmlStats snp.summary.html -csvStats snpEff.snp.csv  -c  ${ref}/ref.config somatic.filter.vcf.gz > somatic.filter.eff.vcf
        bgzip somatic.filter.eff.vcf 
        tabix somatic.filter.eff.vcf.gz
    """
}

    process bamPEforManta {
        publishDir "${params.outdir}/04.bamPE/" , pattern: "*"
        tag "bamPEforManta"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            tuple type,file(bam),ref from bamlist2
        output:
            file "${bam.simpleName}.PE.cram"  into PEbam
            file "${bam.simpleName}.PE.cram.crai" into PEindex
        script:
        """
            samtools view -b --output-fmt CRAM -f 2 ${bam} >${bam.simpleName}.PE.cram
            samtools index ${bam.simpleName}.PE.cram
        """
    }
    process runManta {
        publishDir "${params.outdir}/05.Manta/" , pattern: "*"
        tag "runManta"
        queue "DNA"
        cpus 16
        executor "slurm"
        input:
            file(bam) from PEbam.collect()
            file(bai) from PEindex.collect()
            file normal from normal
            file tumor from tumor
            file(ref) from ref_ch4
        output:
            file "*"
        script:
        if(params.exome){
        """
            configManta.py --referenceFasta=${ref}/ref.fa --bam=${normal.simpleName}.PE.cram --tumorBam=${tumor.simpleName}.PE.cram --runDir Manta --exome
            Manta/runWorkflow.py -j 16 -m local 
        """
        }else{
        """
            configManta.py --referenceFasta=${ref}/ref.fa --bam=${normal.simpleName}.PE.cram --tumorBam=${tumor.simpleName}.PE.cram --runDir Manta --exome
            Manta/runWorkflow.py -j 16 -m local 
        """
        }
    }

    process cram2bam {
        publishDir "${params.outdir}/06.carm2bam/" , pattern: "*"
        tag "bamPEforManta"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            tuple type,file(bam),ref from bamlist3
        output:
            file "${bam.simpleName}.bam"  into bam
            file "${bam.simpleName}.bam.bai" into index
        script:
        """
            samtools view -b --output-fmt BAM  ${bam} >${bam.simpleName}.bam
            samtools index ${bam.simpleName}.bam
        """
    }    
    
    process freec {
        publishDir "${params.outdir}/07.cnv/" , pattern: "*"
        tag "bamPEforManta"
        queue "DNA"
        cpus 8
        executor "slurm"
        input:
            file bam from bam.collect()
            file index from index.collect()
            file normal from normal
            file tumor from tumor
            file(ref) from ref_ch5
        output:
            file "*"
        script:
        if(params.exome){
        """
        mkdir ref
        less -S ${ref}/ref.fa|perl -ne 'if(/>/){@a=split(/\s+/,\$_);\$a[0]=~s/>//g;open Out,">ref/\$a[0].fa"}print Out \$_;'
        echo "[general]"  >? config.txt
        echo "chrLenFile = 03.reference/ref.fa.fai" >> config.txt
        echo "window = 0" >> config.txt
        echo "ploidy = 2" >> config.txt
        echo "breakPointType=4" >> config.txt
        echo "chrFiles =  ref/" >> config.txt
        echo "maxThreads=8" >> config.txt
        echo "breakPointThreshold=1.2" >> config.txt
        echo "readCountThreshold=50" >> config.txt
        echo "[sample]" >> config.txt
        echo "inputFormat = BAM" >>config.txt" >> config.txt
        echo "mateOrientation = RF" >>config.txt" >> config.txt
        echo "[control]" >> config.txt
        echo "inputFormat = BAM" >>config.txt" >> config.txt
        echo "mateOrientation = RF" >>config.txt" >> config.txt
        echo "[target]" >> config.txt
        freec -sample ${tumor.simpleName}.bam -control ${normal.simpleName}.bam -conf config.txt
        Rscript ${baseDir}/bin/significance.R --cnv ${tumor.simpleName}.bam_CNVs --ratio ${tumor.simpleName}.bam_ratio.txt --out ${tumor.simpleName}.cnv_result.txt
        Rscript ${baseDir}/bin//makeGraph.R --ploidy 2 --ratio ${tumor.simpleName}.bam_ratio.txt --out ${tumor.simpleName}
        """
        }else{
        """
        mkdir ref
        less -S ${ref}/ref.fa|perl -ne 'if(/>/){@a=split(/\s+/,\$_);\$a[0]=~s/>//g;open Out,">ref/\$a[0].fa"}print Out \$_;'
        echo "[general]" >config.txt
        echo "chrLenFile = 03.reference/ref.fa.fai" >> config.txt
        echo "ploidy = 2" >> config.txt
        echo "breakPointThreshold = .8" >> config.txt
        echo "window = 50000" >>config.txt
        echo "chrFiles = ref/" >>config.txt
        echo "maxThreads=8" >>config.txt
        echo "[sample]" >>config.txt
        echo "inputFormat = BAM" >>config.txt
        echo "mateOrientation = RF" >>config.txt
        echo "[control]" >> config.txt
        echo "inputFormat = BAM" >> config.txt
        echo "mateOrientation = RF" >> config.txt
        echo "[target]" >>config.txt
        freec -sample ${tumor.simpleName}.bam -control ${normal.simpleName}.bam -conf config.txt
        Rscript ${baseDir}/bin/significance.R --cnv ${tumor.simpleName}.bam_CNVs --ratio ${tumor.simpleName}.bam_ratio.txt --out ${tumor.simpleName}
        Rscript ${baseDir}/bin//makeGraph.R --ploidy 2 --ratio ${tumor.simpleName}.bam_ratio.txt --out ${tumor.simpleName}
        """
        }
    }