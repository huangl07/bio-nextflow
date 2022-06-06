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
library(rrBLUP)
library(dplyr)
library(data.table)
library(pbapply)
setwd(opt$output)
myG<-fread(opt$vcf,head=T,sep="\t")
myY<-read.table(opt$trait,head=T,sep="\t");

transHapmap2numeric=function (G) 
{
    rownamesG <- rownames(G)
    homo <- c("AA", "TT", "CC", "GG")
    hete <- c("AC", "AG", "AT", "CG", "TC", "TG", "CA", "GA", 
        "TA", "GC", "CT", "GT")
    res <- pbapply(G, 2, function(x) {
        SNP<-rep(0,nrow(G))
		SNP[x == "--"]=NA
		SNP[x %in% hete] = 1
		id=names(table(myG$V1));
		id=id[id %in% homo]
		if(length(id) == 1){
			SNP[x == id[1]]=0
		}else{
			SNP[x==id[1]]=0;
			SNP[x==id[2]]=2;
		}
        SNP
    })
    res
    rownames(res) <- rownamesG
    res
}



myG=myG[,-(1:11)]
colname=colnames(myG)
myG=transpose(myG)
myG=as.data.frame(myG)
row.names(myG)=colname
myG$sampleID=row.names(myG)
myG=myG %>% filter(myG$sampleID %in% myY$sampleID) %>% arrange(sampleID)
myY=myY %>% arrange(sampleID)
myG=myG %>% select(!sampleID)
G=transHapmap2numeric(myG)
G=A.mat(G,impute.method="mean",return.imputed=T)
model=mixed.solve(myY[,2],Z=G$imputed)
fix=kin.blup(data=myY,pheno=colnames(myY)[2],K=G$A,geno="sampleID",PEV=T);

pred=as.matrix(G$imputed) %*% as.matrix(model$u)
pred = pred[,1]+ model$beta
result=data.frame(myY,BLUP=fix$g,PEV=fix$PEV,Pred=pred,VarG=fix$Vg,VarE=fix$Ve)

bootstrap<-function(n){
    nrow=nrow(G$imputed)
    subset=sample(1:nrow,ceiling(nrow*0.8))
    G1=G$imput[subset,]
    T1=myY[subset,2]
    pred= as.matrix(G1) %*% as.matrix(model$u )
    pred=pred[,1]+ model$beta
    cor=cor(pred,T1,use="complete")
    df=data.frame(n,cor=cor)
    df
}
cor=pblapply(1:500,bootstrap)
cor=do.call(rbind,cor)
write.table(file=paste(colnames(myY)[2],"fix.cor",sep="."),cor)
write.table(file=paste(colnames(myY)[2],"rrblup.result",sep="."),cor)
save(model,file="model.RDS")

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)