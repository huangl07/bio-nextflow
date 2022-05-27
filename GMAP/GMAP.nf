#!/usr/bin/env nextflow
params.out = "demo"
params.popt="F2"
params.segment=0.5
params.missing=0.3
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
    --chr_only
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
        file "split/compare.sh" into compare_worksh,compare
        path "split"
        file "*"
    script:

    if(params.popt == "CP"){
      if(params.chr){
            """
            perl ${baseDir}/bin/markersplitCP.pl -input ${marker} -fOut scaf.list -dOut split -chr ${chrfile}
            """
        }else{
            """
            perl ${baseDir}/bin/markersplitCP.pl -input ${marker} -fOut scaf.list -dOut split 
            """
        }
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
        if(params.popt == "CP"){
            """
            ln -s ${scaf[0]} ${sca}.bin.marker
            """
        }else{
            if(file("${scaf[0]}").countLines() > 2){
            """
            snpbinner crosspoints -i ${scaf[0]} -o ${sca}.cross -r 0.002
            snpbinner bins -i ${sca}.cross -o ${sca}.bins -l 1000
            perl ${baseDir}/bin/convert2MSTmap.pl -input ${sca}.bins -output ${sca}.bin.marker -popt F2 --chr ${sca} --marker ${scaf[0]}
            """
            }else{
            """
            less -S ${scaf[0]}|sed 's/\t/-/'|sed 's/h/X/g' > ${sca}.bin.marker
            """
            }
        }
}

onlychr=0
if(compare.splitText(by:1).toList().size().val < 1){
   println "haha"
   onlychr=1
   process group_by_chr{
        publishDir "${params.out}/06.premapping", pattern:"*"
        queue "DNA"
        executor "slurm"
        input:
            file bins from binfile1.collect()
        output:
            file "*.lg" into grouping_file
            file "*"
        script:
        """
        cat *.marker > total.markers
        perl ${baseDir}/bin/linkage_by_ref.pl  -i total.markers -o  total.lg    
        """
    }

}else{
    compare_worksh.splitText(by:1).set{compare_para}
    println "test"
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
        """
        perl  ${baseDir}/bin/calculateMLOD.pl  -popt F2 ${para} 
        """
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
                perl ${baseDir}/bin/linkage_by_ref-mlod.pl  -i linkage.mlod.csv -o lg.lg -t 5 -c ${chrfile}
                """
            }else{
                """
                cat *.mlod > Total.mLOD
                python3  ${baseDir}/bin/count_mlod.py -input Total.mLOD -output linkage.mlod.csv 
                perl ${baseDir}/bin/linkage_by_mlod.pl  -i linkage.mlod.csv -k lg -d ./ -n ${params.nchro} -minGroup 1 -b 3 -e 20 
                """
            }
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
                cat *.marker > total.markers
                perl ${baseDir}/bin/splitbyLG-CP.pl  -l ${lg} -i total.markers -d linkagegroups/
            """
        }else{
            if(onlychr > 0){
            """
            cat *.marker > total.markers
            perl ${baseDir}/bin/splitbyLG-NOCP.pl  -l ${lg} -i total.markers -d linkagegroups/ --chr 1
            """
            }else{
            """
            cat *.marker > total.markers
            perl ${baseDir}/bin/splitbyLG-NOCP.pl  -l ${lg} -i total.markers -d linkagegroups/ 
            """
            }
        }
    }


lglist.splitCsv(header:false,sep:'\t').groupTuple().set{lg_list}

process mapping{
    publishDir "${params.out}/07.mapping", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
        tuple lg,marker from lg_list
    output:
        file "*.result.map" into map_result
        file "*.result.csv" into csv_result
        file "*"
    script:
    if(params.popt == "CP"){
        if(params.chr){
            """      
            crosslink_group --inp=${marker[0]} --outbase=${lg}. --knn=30 --matpat_lod=20 --matpat_weights=01P07 --min_lod=6 --redundancy_lod=20    
            ln -s  `wc -l ${lg}.*.loc|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.draw.loc
            cut -f 1 -d " " ${lg}.draw.loc|perl -ne 'chomp;@a=split(/-/,\$_);print \$_,"\\t",\$a[1],"\\n";' > ref.map
            perl ${baseDir}/bin/smooth-CP.pl -l ${lg}.draw.loc -k ${lg} -d ./ -m ref.map   -ami 0.8
            crosslink_group --inp=${lg}.correct.loc --outbase=${lg}.correct. --knn=45 --map=${lg}.correct.
            ln -s  `wc -l ${lg}.correct.*.loc|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.result.csv
            ln -s  `wc -l ${lg}.correct.*.map|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.result.map            
            """
        }else{
            """
            crosslink_group --inp=${marker[0]} --outbase=${lg}. --knn=30 --matpat_lod=20 --matpat_weights=01P07 --min_lod=6 --redundancy_lod=20
            ln -s `wc -l ${lg}.*.loc|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.draw.loc
            crosslink_map --inp=${lg}.draw.loc --map=${lg}.primary.map --out=${lg}.result.csv --ga_gibbs_cycles=10 --ga_iters=300000  --ga_max_hop=1.0 --ga_max_mvseg=1.0 --ga_max_mvdist=1.0 --ga_max_seg=1.0 --gibbs_samples=500  --gibbs_burnin=20 --gibbs_min_prob_1=0.1 --gibbs_min_prob_2=1
            perl ${baseDir}/bin/smooth-CP.pl -l ${lg}.result.csv -k ${lg} -d ./ -m ${lg}.primary.map -win 30 -ami 0.8
            crosslink_group --inp=${lg}.correct.loc --outbase=${lg}.correct. 
            crosslink_map --inp=${lg}.correct.000.loc --map=${lg}.correct.000.map --ga_gibbs_cycles=10 --ga_iters=300000  --ga_max_hop=1.0 --ga_max_mvseg=1.0 --ga_max_mvdist=1.0 --ga_max_seg=1.0 --gibbs_samples=500  --gibbs_burnin=20 --gibbs_min_prob_1=0.1 --gibbs_min_prob_2=1
            ln -s  `wc -l ${lg}.correct.*.loc|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.result.csv
            ln -s  `wc -l ${lg}.correct.*.map|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.result.map 
            """
        }
    }else{
        """
         Rscript ${baseDir}/bin/Asmap.R  --binfile ${marker[0]} --output ${lg}.result.map --popt ${params.popt} --lg ${lg}
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
            cat *.result.csv > total.result.csvs
            perl ${baseDir}/bin/marker2onemap.pl -i total.result.csvs -o total.result.onemap
            perl ${baseDir}/bin/marker2phase.pl -i total.result.csvs -o total.result.phase
            perl ${baseDir}/bin/map-gather.pl -i ./ -o ./
            perl ${baseDir}/bin/mapEstimate.pl -i total.sexAver.map -o total.mapstat
            perl ${baseDir}/bin/drawAligmentRalationMap.pl -m total.sexAver.map -o ./ -k total.phy
            perl ${baseDir}/bin/markerinfo.pl -map total.sexAver.map -input total.result.csvs --pop ${params.popt} --out total
            Rscript ${baseDir}/bin/plotmaps.R --mark total.sexAver.map --out total.map
	        Rscript ${baseDir}/bin/drawbinCP-sexAver.R --mark total.sexAver.phase  --out total.sexAver.bin;

        """
    }else{
        """
         cat *.map > total.maps
         cat *.marker|sort -r|uniq > total.markers
         less -S chr01.result.csv|grep Genotype > tmp.head
         cat *.csv|sort|uniq|grep -v Genotype > tmp.csvs
         cat tmp.head tmp.csvs > total.csvs
         rm tmp.*
         perl ${baseDir}/bin/mapEstimate.pl -i total.maps -o total.mapstat
        perl ${baseDir}/bin/markerinfo.pl -map total.maps -input total.markers --pop ${params.popt} --out total
        Rscript ${baseDir}/bin/drawmap.R --mark total.csvs  --out ./ --pop ${params.popt}
        Rscript ${baseDir}/bin/plotmaps.R --mark total.csvs --out total.map
        Rscript ${baseDir}/bin/drawbinNOCP.R --mark total.csvs  --out ./ --name total
        perl ${baseDir}/bin/drawAligmentRalationMap.pl -m total.maps -o ./ -k total.phy
        """
    }
}

