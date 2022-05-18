#! /mnt/ilustre/app/pub/R/bin/Rscript
library(getopt)
options(bitmapType='cairo')
opt = getopt(matrix(c(
'trt','m',1,'character',
'out','o',2,'character',
'help','h',0,'logical'
),byrow=TRUE, ncol=4));
usage<-function(){
cat("This script is used to plot genetic map
Usage	Rscript map-plot.R for fig 3-9[options]
Options:
	-m, --trt	trt file, forced 
	-o, --out 	out file, forced
	-h, --help	print display this help and exit
")
q(status=1);
}
times<-Sys.time()

if (!is.null(opt$help) ) { usage() }
if (is.null(opt$trt) ) { usage() }
library(dplyr)
library(reshape2)
txt<-read.table(opt$trt,sep="\t",head=T)
id=colnames(txt)[1];
df=melt(txt)
stat=df %>% group_by(variable) %>% summarise(nlevels=nlevels(as.factor(value)))
qname=stat$variable[stat$nlevels > 2]
bname=stat$variable[stat$nlevels == 2]

print(qname)
print(bname)
write.table(txt[,c(id,as.character(qname))],file=paste(opt$out,"qtl.txt",sep="."),quote=FALSE,sep=" ",row.name=F)
if(len(bname) > 0){
write.table(txt[,c(id,as.character(bname))],file=paste(opt$out,"btl.txt",sep="."),quote=FALSE,sep=" ",row.name=F)
}

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
