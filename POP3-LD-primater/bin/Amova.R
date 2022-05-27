#!/usr/bin/env Rscript

times<-Sys.time()

library('getopt')
spec<-matrix(c(
    'help','h',0,"logical",
    'vcf','b',1,"character",
    'group','g',1,'character',
    'out','o',1,'character'
),byrow=TRUE,ncol=4);
opt=getopt(spec)
print_usage<-function(spec=NULL){
    cat(getopt(spec,usage=TRUE));
    cat("Usage example :\n")
    cat("
Usage example:
Usage:
    --vcf input vcf file
    --group input group file
    --out output dir

    --help        usage
\n")
    q(status=1);
}
if(!is.null(opt$help)){print_usage(spec)}
if(is.null(opt$vcf)){print_usage(spec)}
if(is.null(opt$group)){print_usage(spec)}
if(is.null(opt$out)){print_usage(spec)}


if(!dir.exists(opt$out)){dir.create(opt$out)}
library(poppr)
library(pegas)
library("vcfR")
info <- VCFloci(opt$vcf)
popdata<-read.table(opt$group,head=F)
colnames(popdata)=c("sampleID","group");
x <- read.vcf(opt$vcf, from = 1, to = nrow(info))
g<-loci2genind(x)
g<-as.genclone(g)
ploidy(g)<-2
strata(g)<-data.frame(popdata)
g=setPop(g,~group)
setwd(opt$out)

amovacc<-poppr.amova(g,~group)
 write.table(amovacc$componentsofcovariance, sep = ",", quote=F,file = "amova.component.csv")
  write.table(amovacc$result, sep = ",", quote=F,file = "amova.result.csv")

 amovasig   <- randtest(amovacc, nrepet = 999)
 pdf("pop.signif.pdf")
 plot(amovasig)
 dev.off()
 png("pop.signif.png")
 plot(amovasig)
 dev.off()
df=data.frame(names=amovasig$names,obs=amovasig$obs,pvalue=amovasig$pvalue,adjust.p=amovasig$adj.pvalue,alter=amovasig$alter,expvar=amovasig$expvar)
write.table(df,file="amovasig.csv",sep=",",quote=F,row.names=F)
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)