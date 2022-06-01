library('getopt');
times<-Sys.time();

spec = matrix(c(
	'vcf','v',0,'character',
    'pvalue','p',0,'character',
	'trait','t',0,'character',
	'output','o',0,'character',
	'help','c',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript indel_len.r --i  --o  
	
Usage:
	--pop	the input hapmap file
	--trait	the trait file 
	--output	the output dir
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$trait)) { print_usage(spec)}
if ( is.null(opt$vcf)){ print_usage(spec) }
if ( is.null(opt$output)){ print_usage(spec) }
if(!dir.exists(opt$output)){dir.create(opt$output)}
library(GOVS)
library(rrBLUP)
library(dplyr)
myG<-read.table(opt$vcf,head=T,sep="\t")
myY<-read.table(opt$trait,head=T,sep="\t");

myG=myG[,-(1:11)]
myG=t(myG)
myG=as.data.frame(myG)
myG$sampleID=row.names(myG)
myG=myG %>% filter(myG$sampleID %in% myY$sampleID) %>% arrange(sampleID)
myY=myY %>% arrange(sampleID)
myG=myG %>% select(!sampleID)
G=transHapmap2numeric(myG)
G=A.mat(G,impute.method="mean",return.imputed=T)
model=mixed.solve(myY[,2],Z=G$imputed)
fix=kin.blup(data=myY,pheno=colnames(myY)[2],K=G$A,geno="sampleID",PEV=T);

pred=as.matrix(model$u) %*% G$impute + model$beta
fun<-check(g,y){
    pred=   as.matrix(model$u) %*% g + model$beta
    cor=cor(pred,y,use="complete")
    cor
}
fun<-bootstrap(n){
    nrow=nrow(G$imputed)
    subset=samle(1:nrow,ceiling(nrow*0.8))
    G1=G$imput[subset]
    T1=myY[subset,2]
    df=data.frame(n,cor=chek(g,y))
    df
}
cor=lapply(1:500,bootstrap)
cor=do.call(rbind,cor)
result=data.frame(myY,BLUP=fix$g,PEV=fix$PEV,Pred=pred,VarG=fix$Vg,VarE=fix$e)
write.table(fill=paste(colnames(myY)[2],"fix.cor",sep="."),cor)
write.table(fill=paste(colnames(myY)[2],"rrblup.result",sep="."),cor)
save.RDS(model,fill="model.RDS")

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)