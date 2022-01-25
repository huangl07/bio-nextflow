times<-Sys.time();
library('getopt');
options(bitmapType='cairo');
spec = matrix(c(
	'matrix','a',0,'character',
	'group','b',0,'character',
   'output','c',0,'character',
   'contrast','d',0,'character',
	'help' , 'e', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage:
	--matrix	input count matrix file
   --group input group file
   --contrast input condition file [option]
	--output	output dir
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$matrix) ) { print_usage(spec) }
if ( is.null(opt$group) ) { print_usage(spec) }
if ( is.null(opt$output) ) { print_usage(spec) }
library(DESeq2)
library("BiocParallel")
register(MulticoreParam(8))
library(ggplot2)

group<-read.table(opt$group,head=T,comment.char="^",sep="\t")
gene<-read.csv(opt$matrix,head=T)
if(!is.null(opt$contrast)){condition=read.table(opt$contrast,head=T);}
if(!dir.exists(opt$output)){dir.create(opt$output)}
setwd(opt$output)
rownames(gene)=gene[,1]
colnames(group)=c("sampleID","group")
rownames(group)=group[,1]
dds<-DESeqDataSetFromMatrix(countData=gene[,-1],colData=group,design =~ group)
dds <- DESeq(dds,parallel=T)

 contrast<-function(i){
      tag=paste(condition[i,1],condition[i,2],sep="_vs_")
      print(tag)
      res <- results(dds,contrast=c("group",condition[i,1],condition[i,2]),parallel=T)
      res=as.data.frame(res)
      res$padj[is.na(res$padj)]=1
      res$type="nosig";
      res$type[res$log2FoldChange > 1 & res$padj < 0.05]="up";
      res$type[res$log2FoldChange < -1 & res$padj < 0.05]="down";
      write.table(file=paste(tag,".xls",sep=""),res)
      write.table(file=paste(tag,"sig","xls",sep=.""),res[res$type != "nosig",]);
      p=ggplot(na.omit(res))+geom_point(aes(x=log2FoldChange,y=-1*log2(padj),col=type))
      ggsave(file=paste(tag,".volcano.pdf",sep=""),p,device="pdf")
   }
if (is.null(opt$contrast)){
   lgroup=levels(as.factor(group$group))
   condition=t(combn(levels(as.factor(group$group)),2))
   lapply(c(1:ncol(condition)),contrast)
}else{
   lapply(c(1:nrow(condition)),contrast)
}


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()


