#!/usr/bin/env Rscript

times<-Sys.time()

library('getopt')
spec<-matrix(c(
    'help','h',0,"logical",
    'vcf','b',1,"character",
    'group','g',1,'character',
    'stack','stack',1,'character',
    'pic','p',1,'character',
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
    --stack input stack file
    --out output dir
    --pic   input pic result
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
library(dplyr)
library(PopGenome)
library(hierfstat)
info <- VCFloci(opt$vcf)
popdata<-read.table(opt$group,head=F)
colnames(popdata)=c("sampleID","group");
x <- read.vcf(opt$vcf, from = 1, to = nrow(info))
g<-loci2genind(x)
g<-as.genclone(g)
ploidy(g)<-2
strata(g)<-data.frame(popdata)
g=setPop(g,~group)
diversity <- poppr(g)
print("poppr down")
f=genind2hierfstat(g)
gid=levels(f$pop)
stat<-function(i){df=basic.stats(f[f$pop==gid[i]]);r=df$result;df1=data.frame(t(df$overall),Pop=gid[i]);df1}
dfs=lapply(1:length(gid),stat)
dfs1=do.call(rbind,dfs)
diversity=left_join(diversity,dfs1,by="Pop")
print("hierfstat down")


VCF_split_into_scaffolds(opt$vcf,"scaffolds")
popGenome=readData("scaffolds/",format="VCF")
list=split(popdata$sampleID,f=as.factor(popdata$group))
pop=set.populations(popGenome,list)
neutrality=get.neutrality(neutrality.stats(pop))
merge<-function(i){df=data.frame(neutrality[[i]]);df$Pop=names(list)[i];df}
dfs=lapply(1:length(list),merge)
result=do.call(rbind,dfs)
result[is.na(result)]=0
result=result %>% group_by(Pop) %>% summarise(Fu.Li.F=mean(Fu.Li.F),Fu.Li.D=mean(Fu.Li.D),Rozas.R_2=mean(Rozas.R_2),Fu.F_S=mean(Fu.F_S),Fay.Wu.H=mean(Fay.Wu.H),Zeng.E=mean(Zeng.E),Strobeck.S=mean(Strobeck.S))
ndf=left_join(diversity,result,by="Pop")

if(!is.null(opt$pic)){
    pic<-read.table(opt$pic)
    colnames(pic)=c("Pop","PIC")
    ndf=left_join(ndf,pic,by="Pop")
}
df1=read.talbe(opt$stack,sep="\t",head=T)
colnames(df1)[1]="Pop";
left_join(ndf,df1,by="Pop")
ndf$Ne=diversity$pi * 4 /1e-8




write.table(ndf,file=paste(opt$out,"diversity.csv",sep="/"),row.names=F,sep=",")

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)