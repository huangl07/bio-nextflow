#!/usr/bin/env Rscript
times<-Sys.time();
library('getopt');
options(bitmapType='cairo');
spec = matrix(c(
	'cnv','c',0,'character',
	'ratio','r',0,'character',
	'out','o',0,'character',
	'help' , 'e', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	
Usage:
    --cnv input cnv file
    --ratio input ratio file
    --out input out file
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$cnv) ) { print_usage(spec) }
if ( is.null(opt$ratio) ) { print_usage(spec) }
if ( is.null(opt$out) ) { print_usage(spec) }

library(rtracklayer)

dataTable <-read.table(opt$ratio, header=TRUE);
ratio<-data.frame(dataTable)

dataTable <- read.table(opt$cnv, header=FALSE)
cnvs<- data.frame(dataTable)

ratio$Ratio[which(ratio$Ratio==-1)]=NA

cnvs.bed=GRanges(cnvs[,1],IRanges(cnvs[,2],cnvs[,3]))  
ratio.bed=GRanges(ratio$Chromosome,IRanges(ratio$Start,ratio$Start),score=ratio$Ratio)

overlaps <- subsetByOverlaps(ratio.bed,cnvs.bed)
normals <- setdiff(ratio.bed,cnvs.bed)
normals <- subsetByOverlaps(ratio.bed,normals)

#mu <- mean(score(normals),na.rm=TRUE)
#sigma<- sd(score(normals),na.rm=TRUE)

#hist(score(normals),n=500,xlim=c(0,2))
#hist(log(score(normals)),n=500,xlim=c(-1,1))

#shapiro.test(score(normals)[which(!is.na(score(normals)))][5001:10000])
#qqnorm (score(normals)[which(!is.na(score(normals)))],ylim=(c(0,10)))
#qqline(score(normals)[which(!is.na(score(normals)))], col = 2)

#shapiro.test(log(score(normals))[which(!is.na(score(normals)))][5001:10000])
#qqnorm (log(score(normals))[which(!is.na(score(normals)))],ylim=(c(-6,10)))
#qqline(log(score(normals))[which(!is.na(score(normals)))], col = 2)

numberOfCol=length(cnvs)

for (i in c(1:length(cnvs[,1]))) {
  values <- score(subsetByOverlaps(ratio.bed,cnvs.bed[i]))
  #wilcox.test(values,mu=mu)
  W <- function(values,normals){resultw <- try(wilcox.test(values,score(normals)), silent = TRUE)
	if(class(resultw)=="try-error") return(list("statistic"=NA,"parameter"=NA,"p.value"=NA,"null.value"=NA,"alternative"=NA,"method"=NA,"data.name"=NA)) else resultw}
  KS <- function(values,normals){resultks <- try(ks.test(values,score(normals)), silent = TRUE)
	if(class(resultks)=="try-error") return(list("statistic"=NA,"p.value"=NA,"alternative"=NA,"method"=NA,"data.name"=NA)) else resultks}
  #resultks <- try(KS <- ks.test(values,score(normals)), silent = TRUE)
  #	if(class(resultks)=="try-error") NA) else resultks
  cnvs[i,numberOfCol+1]=W(values,normals)$p.value
  cnvs[i,numberOfCol+2]=KS(values,normals)$p.value
  }

if (numberOfCol==5) {
  names(cnvs)=c("chr","start","end","copy number","status","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue")  
}
if (numberOfCol==7) {
  names(cnvs)=c("chr","start","end","copy number","status","genotype","uncertainty","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue")  
}
if (numberOfCol==9) {
  names(cnvs)=c("chr","start","end","copy number","status","genotype","uncertainty","somatic/germline","precentageOfGermline","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue")  
}
write.table(cnvs, file=paste(opt$out,"cnv.result",sep="."),sep="\t",quote=F,row.names=F)