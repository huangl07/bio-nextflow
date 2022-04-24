#!/usr/bin/env nextflow
params.out = "demo"
params.popt="F2"
params.segment=0.05
params.missing=0.5
params.Pdep=10
params.Odep=2
params.nchro=1
params.help = false
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
   

    --vcf   <file>  input vcf file
    --popt  <str>   population type CP/BCi/Fi/Rix/
    --out   <dir>   output dir
    --p1    <str>   input p1 IDs
    --p2    <str>   input p2 IDs
    --Pdep  <num>   input parents depth
    --Odep  <num>   input offspring depth
    --segment   <num>   segment 0.05
    --missing   <num>   missing 0.5
    --chr   <file>  chr.list
    --nchro <num>   chr number
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}

vcf=file(params.vcf)
chr=file(params.chr)
process vcf2table{
    publishDir "${params.out}/01.vcfconvert", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file vcf from vcf
    output:
        file "pop.filtered.marker" into filtered
        file "*"
    script:
    """
        perl ${baseDir}/bin/vcf2marker.pl -vcf ${vcf} -out pop.primary.marker -PID ${params.p1} -MID ${params.p2} -Pdep ${params.Pdep} -Odep ${params.Odep} 
        perl ${baseDir}/bin/markerfilter.pl -input pop.primary.marker -output pop  -mis ${params.missing} -seg ${params.segment} -popt ${params.popt}
    """
}

process markersplit{
    publishDir "${params.out}/02.markersplit", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file marker from filtered
        file chrfile from chr
    output:
        file "scaf.list" into scaflist
        file "split/compare.sh" into compare_worksh
        path "split"
        file "*"
    script:

    if(params.popt == "CP"){
    """

    """
    }else{
        if(params.chr){
            """
            perl ${baseDir}/bin/markersplitNOCP.pl -input ${marker} -fOut scaf.list -dOut split -chr ${chrfile}
            """
        }else{
            """
            perl ${baseDir}/bin/markersplitNOCP.pl -input ${marker} -fOut scaf.list -dOut split 
            """
        }
    }
}


scaflist.splitCsv(header:false,sep:'\t').groupTuple().set{scaf_list}

process binmap{
    publishDir "${params.out}/03.bin", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        tuple val(sca),scaf from scaf_list
    output:
        file "*.bin.marker" into binfile1,binfile2
        file "*"
    script:
    if(file("${scaf[0]}").countLines() > 2){
        if(params.popt == "CP"){
        """

        """
        }else{
        """

        snpbinner crosspoints -i ${scaf[0]} -o ${sca}.cross -r 0.002
        snpbinner bins -i ${sca}.cross -o ${sca}.bins -l 1000
        perl ${baseDir}/bin/convert2MSTmap.pl -input ${sca}.bins -output ${sca}.bin.marker -popt F2 --chr ${sca} --marker ${scaf[0]}

        """
        }
    }else{
        """
        less -S ${scaf[0]}|sed 's/\t/-/'|sed 's/h/X/g' > ${sca}.bin.marker
        """
    }
}
compare_worksh.splitText(by:1).set{compare_para}

process mlodcalc{
    publishDir "${params.out}/04.calc", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        val para from compare_para
        file binfile from binfile1.collect()
    output:
        file "*.mlod" into mlod_file
        file "*"
    script:

    if(params.popt == "CP"){
    """

    """
    }else{
    """
        perl  ${baseDir}/bin/calculateMLOD.pl  -popt F2 ${para} 
    """
    }
}


process grouping{
    publishDir "${params.out}/05.group", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file mlod from mlod_file.collect()
        file chrfile from chr
    output:
        file "*.lg" into grouping_file
        file "*"
    script:
    if(params.chr){
        """
        cat *.mlod > Total.mLOD
        python3  ${baseDir}/bin/count_mlod.py -input Total.mLOD -output linkage.mlod.csv
        perl ${baseDir}/bin/linkage_by_ref.pl  -i linkage.mlod.csv -o lg.lg -t 5 -c ${chrfile}
        """
    }else{
        """
        cat *.mlod > Total.mLOD
        python3  ${baseDir}/bin/count_mlod.py -input Total.mLOD -output linkage.mlod.csv 
        perl ${baseDir}/bin/linkage_by_mlod.pl  -i linkage.mlod.csv -k lg -d ./ -n ${params.nchro} -minGroup 1 -b 3 -e 20 
        """
    }
}


process premapping{
    publishDir "${params.out}/06.premapping", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file lg from grouping_file
        file bins from binfile2.collect()
    output:
        file "linkagegroups/lg.list" into lglist
        file "linkagegroups/*.marker" into lgmarkers1,lgmakrers2
        file "*"
    script:
    if(params.popt == "CP"){
        """
      
        """
    }else{
        """
        cat *.marker > total.markers
        perl ${baseDir}/bin/splitbyLG-NOCP.pl  -l ${lg} -i total.markers -d linkagegroups/
        """
    }
}
lglist.splitCsv(header:false,sep:'\t').groupTuple().set{lg_list}

process mapping{
    publishDir "${params.out}/07.mapping", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        tuple lg,marker from lg_list
        file lg_markers from lgmarkers1
    output:
        file "${lg}.map" into map_result
        file "${lg}.csv" into csv_result
        file "*"
    script:
    if(params.popt == "CP"){
        """
      
        """
    }else{
        """
         Rscript ${baseDir}/bin/Asmap.R  --binfile ${marker[0]} --output ${lg}.map --popt ${params.popt} --lg ${lg}
        """
    }
}

process mapEvaluate{
    publishDir "${params.out}/08.evaluate", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        file maps from map_result.collect()
        file lg from lgmakrers2.collect()
        file csv from csv_result.collect()
    output:
        file "*"
    script:
    if(params.popt == "CP"){
        """
      
        """
    }else{
        """
         cat *.map > total.maps
         cat *.marker|sort -r|uniq > total.markers
         cat *.csv|sort -r|uniq > total.csvs
        perl ${baseDir}/bin/markerinfo.pl -map total.maps -input total.markers --pop ${params.popt} --out total
        Rscript ${baseDir}/bin/drawmap.R --mark total.csvs  --out ./ --pop ${params.popt}
        Rscript ${baseDir}/bin/plotmaps.R --mark total.csvs --out total.map
        Rscript ${baseDir}/bin/drawbinNOCP.R --mark total.csvs  --out ./ --name total
        perl ${baseDir}/bin/drawAligmentRalationMap.pl -m total.maps -o ./ -k total.phy
        """
    }
}