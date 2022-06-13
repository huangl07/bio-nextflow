#!/usr/bin/env nextflow
params.fqlist = 'fq.list'
params.outdir = "demo"
params.help = false
params.method="data2"
params.meta="meta.list"
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --fqlist    <file>  rawdata file list
    --outdir    <dir>   output dir
    --method    <str> data2 or debular
    --meta  <file> group list
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

fq_list = file(params.fqlist)
group_list=file(params.meta)

process importfq{
    publishDir "${params.outdir}/01.importfq", pattern:"*"
    queue "DNA"
    executor "slurm"
    beforeScript 'source /mnt/ilustre/users/dna/.bash_conda'
    conda '/mnt/ilustre/users/dna/miniconda3/envs/qiime2-2021.8'
    input:
        file fqlist from fq_list
    output:
        file "*"
        file "paired-end-demux.qza" into importdata1,importdata2
    script:
    
    """
        qiime tools import \
            --type \'SampleData[PairedEndSequencesWithQuality]\' \
            --input-path ${fqlist} \
            --output-path paired-end-demux.qza \
            --input-format PairedEndFastqManifestPhred33V2
    """
}

process DeNoise{
    publishDir "${params.outdir}/02.DeNoise", pattern:"*"
    queue "DNA"
    executor "slurm"
    if(params.method == "data2"){
    cpus 8
    }else{
    cpus 1
    }
    beforeScript 'source /mnt/ilustre/users/dna/.bash_conda'
    conda '/mnt/ilustre/users/dna/miniconda3/envs/qiime2-2021.8'
    input:
        file demux from importdata1
    output:
        file "*"
        file "table.qza" into tableqza1,tableqza2,table
        file "rep-seqs.qza" into repseqs1,repseqs2,seq
        //file "paired-end-demux.qza" into importdata
    script:
    if (params.method == "data2"){
        """
            qiime dada2 denoise-single \
            --i-demultiplexed-seqs ${demux} \
            --p-trim-left 0 \
            --p-trunc-len 120 \
            --p-n-threads 8 \
            --o-representative-sequences rep-seqs-dada2.qza \
            --o-table table-dada2.qza \
            --o-denoising-stats stats-dada2.qza
            qiime tools export --input-path  stats-dada2.qza --output-path ./
            mv rep-seqs-dada2.qza rep-seqs.qza
            mv table-dada2.qza table.qza
        """
    }else{
         """
        qiime quality-filter q-score \
            --i-demux demux.qza \
            --o-filtered-sequences demux-filtered.qza \
            --o-filter-stats demux-filter-stats.qza
        qiime deblur denoise-16S \
            --i-demultiplexed-seqs demux-filtered.qza \
            --p-trim-length 120 \
            --o-representative-sequences rep-seqs-deblur.qza \
            --o-table table-deblur.qza \
            --p-sample-stats \
            --o-stats deblur-stats.qza\
        mv rep-seqs-deblur.qza rep-seqs.qza
        mv table-deblur.qza table.qza
        """
    }
}


process classify{
    publishDir "${params.outdir}/03.Classify", pattern:"*"
    queue "DNA"
    executor "slurm"
    cpus 8
    beforeScript 'source /mnt/ilustre/users/dna/.bash_conda'
    conda '/mnt/ilustre/users/dna/miniconda3/envs/qiime2-2021.8'
    input:
        file reads from repseqs1
        file table from tableqza1
        file group from group_list
    output:
        file "*"
        file "taxonomy.qza" into taxonomy
    script:
    
    """
        qiime feature-classifier classify-sklearn \
          --i-classifier /mnt/ilustre/users/dna/database/silva-138-99-nb-classifier.qza\
          --i-reads ${reads} \
          --o-classification taxonomy.qza \
          --p-n-jobs 8
        qiime taxa barplot \
            --i-table ${table} \
            --i-taxonomy taxonomy.qza \
            --m-metadata-file ${group} \
            --o-visualization taxa-bar-plots.qzv
        qiime metadata tabulate \
            --m-input-file taxonomy.qza \
            --o-visualization taxonomy.qzv
        qiime tools export --input-path taxonomy.qzv --output-path ./
    """

}

process Tree{
    publishDir "${params.outdir}/04.Tree", pattern:"*"
    queue "DNA"
    executor "slurm"
    beforeScript 'source /mnt/ilustre/users/dna/.bash_conda'
    conda '/mnt/ilustre/users/dna/miniconda3/envs/qiime2-2021.8'
    input:
        file repseq from repseqs2
    output:
        file "*"
        file("rooted-tree.qza") into rooted_tree
        file("unrooted_tree.nwk") into unrooted_tree
    script:
    
    """
        qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences ${repseq} \
        --o-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza \
        --o-tree unrooted-tree.qza \
        --o-rooted-tree rooted-tree.qza
        qiime tools export --input-path rooted-tree.qza --output-path  ./
        mv tree.nwk rooted_tree.nwk
        qiime tools export --input-path unrooted-tree.qza --output-path  ./
        mv tree.nwk unrooted_tree.nwk
    """
}


process Rarecurve{
    publishDir "${params.outdir}/06.Rarecurve", pattern:"*"
    queue "DNA"
    executor "local"
    input:
        file table from table
        file rooted from rooted_tree
        file sequences from seq
        file taxon from taxonomy
        file group from group_list
    output:
        file "*"
    script:
    
    """
    #!/usr/bin/env Rscript
    library(MicrobiotaProcess)
    library(phyloseq)
    library(tidyverse)
    library(RColorBrewer)
    otu <- "${table}"
    rep <- "${sequences}"
    tree <- "${rooted}"
    tax <- "${taxon}"
    sample <- "${group}"
    ps <- mp_import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,mapfilename=sample,treeqza=tree)
    ps %<>% mp_rrarefy()
    ps1=ps %>%  mp_cal_rarecurve(
        .abundance = RareAbundance,
        chunks = 400
    )
    p1 <- ps1 %>% 
      mp_plot_rarecurve(
        .rare = RareAbundanceRarecurve, 
        .alpha = c("Observe","Chao1"),
      )
    ggsave(file="RareAbundanceRarecurve.pdf",device="pdf",p1)
    """
}
process Alpha{
    publishDir "${params.outdir}/07.Alpha", pattern:"*"
    queue "DNA"
    executor "local"
    input:
        file table from table
        file rooted from rooted_tree
        file sequences from seq
        file taxon from taxonomy
        file group from group_list
    output:
        file "*"
    script:
    
    """
    #!/usr/bin/env Rscript
    library(MicrobiotaProcess)
    library(phyloseq)
    library(tidyverse)
    library(RColorBrewer)
    otu <- "${table}"
    rep <- "${sequences}"
    tree <- "${rooted}"
    tax <- "${taxon}"
    sample <- "${group}"
    ps <- mp_import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,mapfilename=sample,treeqza=tree)
    ps %<>% mp_rrarefy()
    ps1=ps %>%  mp_cal_alpha(
        .abundance = RareAbundance,
        chunks = 400
    )
    write.table(file="alpha.index",ps1)

    """
}
process Abundance{
    publishDir "${params.outdir}/07.Abundance", pattern:"*"
    queue "DNA"
    executor "local"
    input:
        file table from table
        file rooted from rooted_tree
        file sequences from seq
        file taxon from taxonomy
        file group from group_list
    output:
        file "*"
    script:
    
    """
    #!/usr/bin/env Rscript
    library(MicrobiotaProcess)
    library(phyloseq)
    library(tidyverse)
    library(RColorBrewer)
    otu <- "${table}"
    rep <- "${sequences}"
    tree <- "${rooted}"
    tax <- "${taxon}"
    sample <- "${group}"
    ps <- mp_import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,mapfilename=sample,treeqza=tree)
    ps %<>% mp_rrarefy()
    ps1=ps %>%  mp_plot_abundance(
           .abundance=RareAbundance,
           .group=group, 
           taxa.class = Phylum, 
           topn = 20
         )
    ggsave(file="abundance.pdf",ps1)

    """
}
process Dist{
    publishDir "${params.outdir}/08.Dist", pattern:"*"
    queue "DNA"
    executor "local"
    input:
        file table from table
        file rooted from rooted_tree
        file sequences from seq
        file taxon from taxonomy
        file group from group_list
    output:
        file "*"
    script:
    
    """
    #!/usr/bin/env Rscript
    library(MicrobiotaProcess)
    library(phyloseq)
    library(tidyverse)
    library(RColorBrewer)
    otu <- "${table}"
    rep <- "${sequences}"
    tree <- "${rooted}"
    tax <- "${taxon}"
    sample <- "${group}"
    ps <- mp_import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,mapfilename=sample,treeqza=tree)
    ps %<>% mp_rrarefy()
    ps %<>% mp_decostand(.abundance=Abundance)
    ps1=ps %>% mp_cal_dist(.abundance=hellinger,distmethod="bray")
    p2 <- ps1 %>% mp_plot_dist(.distmethod = bray, .group = group)
    ggsave(file="Distance.pdf",p2)

    """
}

process Hierarchical {
    publishDir "${params.outdir}/07.Hierarchical", pattern:"*"
    queue "DNA"
    executor "local"
    input:
        file table from table
        file rooted from rooted_tree
        file sequences from seq
        file taxon from taxonomy
        file group from group_list
    output:
        file "*"
    script:
    
    """
    #!/usr/bin/env Rscript
    library(MicrobiotaProcess)
    library(phyloseq)
    library(tidyverse)
    library(RColorBrewer)
    library(ggtree)
    library(ggtreeExtra)
    otu <- "${table}"
    rep <- "${sequences}"
    tree <- "${rooted}"
    tax <- "${taxon}"
    sample <- "${group}"
    ps <- mp_import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,mapfilename=sample,treeqza=tree)
    ps %<>% mp_rrarefy()
    ps %<>% mp_decostand(.abundance=Abundance)
    ps %<>% mp_cal_abundance(add=on)
    ps %<>% mp_cal_clust(
         .abundance=hellinger, 
         distmethod = "bray",
         hclustmethod = "average", # (UPGAE)
         action = "add" # action is used to control which result will be returned
       )
    sample.clust <- ps %>% mp_extract_internal_attr(name='SampleClust')
    phyla.tb <- ps %>% mp_extract_abundance(taxa.class=Phylum, topn=10)
    phyla.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Phyla="label")
    p <- ggtree(sample.clust) + 
       geom_tippoint(aes(color=group)) +
       geom_tiplab()
    p1 <- p + 
      geom_fruit(
         data=phyla.tb,
         geom=geom_col,
         mapping = aes(x = RelRareAbundanceBySample, 
                       y = Sample, 
                       fill = Phyla,
                 ),
         orientation = "y",
         offset = 0.4,
         pwidth = 3, 
         axis.params = list(axis = "x", 
                            title = "The relative abundance of phyla (%)",
                            title.size = 4,
                            text.size = 2, 
                            vjust = 1),
         grid.params = list()
      )
    ggsave(file="Hierarchical.pdf",p1)

    """
}
process PCOA{
    publishDir "${params.outdir}/07.PCOA", pattern:"*"
    queue "DNA"
    executor "local"
    input:
        file table from table
        file rooted from rooted_tree
        file sequences from seq
        file taxon from taxonomy
        file group from group_list
    output:
        file "*"
    script:
    
    """
    #!/usr/bin/env Rscript
    library(MicrobiotaProcess)
    library(phyloseq)
    library(tidyverse)
    library(RColorBrewer)
    otu <- "${table}"
    rep <- "${sequences}"
    tree <- "${rooted}"
    tax <- "${taxon}"
    sample <- "${group}"
    ps <- mp_import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,mapfilename=sample,treeqza=tree)
    ps %<>% mp_rrarefy()
    ps %<>% mp_decostand(.abundance=Abundance)
    ps %<>% mp_cal_pcoa(.abundance=hellinger, distmethod="bray")
    p1 <- ps %>%
        mp_plot_ord(
          .ord = pcoa, 
          .group = group, 
          .color = group, 
          .size = 1.2,
          .alpha = 1,
          ellipse=TRUE
        )
    ggsave(file="pcoa.pdf",p1)

    """
}
process adonis{
    publishDir "${params.outdir}/07.adnonis", pattern:"*"
    queue "DNA"
    executor "local"
    input:
        file table from table
        file rooted from rooted_tree
        file sequences from seq
        file taxon from taxonomy
        file group from group_list
    output:
        file "*"
    script:
    """
    #!/usr/bin/env Rscript
    library(MicrobiotaProcess)
    library(phyloseq)
    library(tidyverse)
    library(RColorBrewer)
    otu <- "${table}"
    rep <- "${sequences}"
    tree <- "${rooted}"
    tax <- "${taxon}"
    sample <- "${group}"
    ps <- mp_import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,mapfilename=sample,treeqza=tree)
    ps %<>% mp_rrarefy()
     ps %<>% mp_decostand(.abundance=Abundance)
    ps %<>% mp_cal_pcoa(.abundance=hellinger, distmethod="bray")
    ps %<>% mp_adonis(.abundance=hellinger, .formula=~group, distmethod="bray", permutations=9999, action="add")
    df=ps %>% mp_extract_internal_attr(name=adonis)
    write.table(file="adonis.table",df[['aov.tab']])
    """
}
process lafse{
    publishDir "${params.outdir}/07.lafse", pattern:"*"
    queue "DNA"
    executor "local"
    input:
        file table from table
        file rooted from rooted_tree
        file sequences from seq
        file taxon from taxonomy
        file group from group_list
    output:
        file "*"
    script:
    """
    #!/usr/bin/env Rscript
    library(ggtree)
    library(ggtreeExtra)
    library(ggplot2)
    library(MicrobiotaProcess)
    library(tidytree)
    library(ggstar)
    library(forcats)
    otu <- "${table}"
    rep <- "${sequences}"
    tree <- "${rooted}"
    tax <- "${taxon}"
    sample <- "${group}"
    ps <- mp_import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,mapfilename=sample,treeqza=tree)
    ps %<>% mp_rrarefy()
    ps %<>%
    mp_diff_analysis(
       .abundance = RelRareAbundanceBySample,
       .group = group,
       first.test.alpha = 0.01
    )
    taxa.tree <- ps %>% mp_extract_tree(type="taxatree")
    taxa.tree %>% select(label, nodeClass, LDAupper, LDAmean, LDAlower, Sign_time, pvalue, fdr) %>% dplyr::filter(!is.na(fdr))
    p1 <- ggtree(
            taxa.tree,
            layout="radial",
            size = 0.3
        ) +
        geom_point(
            data = td_filter(!isTip),
            fill="white",
            size=1,
            shape=21
        )
    p2 <- p1 +
        geom_hilight(
            #data = td_filter(nodeClass == "Phylum"), # This is supported by ggtree >=3.1.3
            mapping = aes(subset = nodeClass == "Phylum", 
                        node = node, 
                        fill = label)
        )
    p3 <- p2 +
        ggnewscale::new_scale("fill") +
        geom_fruit(
            data = td_unnest(RareAbundanceBySample),
            geom = geom_star,
            mapping = aes(
                        x = fct_reorder(Sample, group, .fun=min),
                        size = RelRareAbundanceBySample,
                        fill = group,
                        subset = RelRareAbundanceBySample > 0
                    ),
            starshape = 13,
            starstroke = 0.25,
            offset = 0.04,
            pwidth = 0.8,
            grid.params = list(linetype=2)
        ) +
        scale_size_continuous(
            name="Relative Abundance (%)",
            range = c(1, 3)
        ) 
    p4 <- p3 + geom_tiplab(size=2, offset=7.2)
    p5 <- p4 +
      ggnewscale::new_scale("fill") +
      geom_fruit(
         geom = geom_col,
         mapping = aes(
                       x = LDAmean,
                       fill = Sign_time,
                       subset = !is.na(LDAmean)
                       ),
         orientation = "y",
         offset = 0.3,
         pwidth = 0.5,
         axis.params = list(axis = "x",
                            title = "Log10(LDA)",
                            title.height = 0.01,
                            title.size = 2,
                            text.size = 1.8,
                            vjust = 1),
         grid.params = list(linetype = 2)
      )
    p6 <- p5 +
        ggnewscale::new_scale("size") +
        geom_point(
            data=td_filter(!is.na(fdr)),
            mapping = aes(size = -log10(fdr),
                        fill = Sign_time,
                        ),
            shape = 21,
        )
    p=p6 + theme(
            legend.key.height = unit(0.3, "cm"),
            legend.key.width = unit(0.3, "cm"),
            legend.spacing.y = unit(0.02, "cm"),
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 9),
            )
    ggsave(file="LEfSe.pdf",p2)
    """
}
process diversity{
    publishDir "${params.outdir}/06.Diversity", pattern:"*"
    queue "DNA"
    executor "slurm"
    beforeScript 'source /mnt/ilustre/users/dna/.bash_conda'
    conda '/mnt/ilustre/users/dna/miniconda3/envs/qiime2-2021.8'
   input:
       file table from tableqza2
       file tree from rooted_tree
       file group from group_list
   output:
       file "*"
   script:
   
   """
   qiime diversity core-metrics-phylogenetic \
       --i-phylogeny ${tree} \
       --i-table ${table}\
       --p-sampling-depth 1103 \
       --output-dir core-metrics-results \
       --m-metadata-file ${group}
   qiime diversity alpha-group-significance \
       --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
       --m-metadata-file ${group} \
       --o-visualization core-metrics-results/faith-pd-group-significance.qzv
   qiime diversity alpha-group-significance \
       --i-alpha-diversity core-metrics-results/evenness_vector.qza \
       --m-metadata-file ${group} \
       --o-visualization core-metrics-results/evenness-group-significance.qzv
   qiime diversity beta-group-significance \
       --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
       --m-metadata-file ${group} \
       --m-metadata-column body-site \
       --o-visualization core-metrics-results/unweighted-unifrac-body-site-significance.qzv \
       --p-pairwise
   qiime diversity beta-group-significance \
       --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
       --m-metadata-file ${group} \
       --m-metadata-column subject \
       --o-visualization core-metrics-results/unweighted-unifrac-subject-group-significance.qzv \
       --p-pairwise
   qiime emperor plot \
       --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
       --m-metadata-file ${group} \
       --p-custom-axes days-since-experiment-start \
       --o-visualization core-metrics-results/unweighted-unifrac-emperor-days-since-experiment-start.qzv
   qiime emperor plot \
       --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
       --m-metadata-file ${group} \
       --p-custom-axes days-since-experiment-start \
       --o-visualization core-metrics-results/bray-curtis-emperor-days-since-experiment-start.qzv
   qiime diversity alpha-rarefaction \
           --i-phylogeny ${tree} \
           --i-table ${table}\
           --p-max-depth 4000 \
           --m-metadata-file sample-metadata.tsv \
           --o-visualization alpha-rarefaction.qzv
   """
}
