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
	--vcf	the input hapmap file
	--trait	the trait file 
	--output	the output dir
    --pvalue
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$trait)) { print_usage(spec)}
if ( is.null(opt$vcf)){ print_usage(spec) }
if ( is.null(opt$output)){ print_usage(spec) }
if(!dir.exists(opt$output)){dir.create(opt$output)}
opt$pvalue=as.numeric(opt$pvalue)
library(dplyr)
library(rMVP)
opt$vcf=normalizePath(opt$vcf)
opt$trait=normalizePath(opt$trait)
setwd(opt$output);
MVP.Data(fileVCF=opt$vcf,filePhe=opt$trait,sep.phe="\t",SNP.eff="Add",out="mvp.hmp",priority="memory",ncpus=8)
genotype <- attach.big.matrix("mvp.hmp.geno.desc")
phenotype <- read.table("mvp.hmp.phe",head=TRUE)
pheid=colnames(phenotype)[2]
map <- read.table("mvp.hmp.geno.map" , head = TRUE)
#MVP.Hist(phe=phenotype[,c(1,i)], file="png", breakNum=30, dpi=300)
imMVP <- MVP(
		phe=phenotype,
		geno=genotype,
		map=map,
		nPC.GLM=3,
		nPC.MLM=3,
		nPC.FarmCPU=3,
		priority="speed",
		ncpus=8,
		vc.method="EMMA",
		maxLoop=10,
		method.bin="FaST-LMM",
		p.threshold=opt$pvalue,
		file.type="pdf",
		method=c("GLM", "MLM", "FarmCPU")
)
lamda<-function(qvalue=NULL){
	z=qnorm(qvalue/2)
	return(round(median(z^2, na.rm = TRUE) / qchisq(0.5, 1), 3))
}
colnames(imMVP$glm.results)[1:2]=c("G.effect","G.se")
colnames(imMVP$mlm.results)[1:2]=c("M.effect","M.se")
colnames(imMVP$farmcpu.results)[1:2]=c("F.effect","F.se")
iMVP.res <- cbind(map,imMVP$glm.results, imMVP$mlm.results, imMVP$farmcpu.results)
imMVP<-select(iMVP.res,1:3|ends_with("MLM") | ends_with("GLM")|ends_with("FarmCPU"))
lamdas<-apply(imMVP[,4:ncol(imMVP)],2,lamda)
lamdas=as.data.frame(lamdas)
iMVP.res$threshold=opt$pvalue/nrow(iMVP.res)
write.table(file=paste(pheid,"lamda.csv",sep="."),lamdas,sep="\t",quote=F,row.names=F)
write.table(file=paste(pheid,"GWAS.result.csv",sep="."),iMVP.res,sep="\t",quote=F,row.names=F)
#MVP.Report(imMVP, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=c(1e-6,1e-4),threshold.lty=c(1,2),col=c("grey60","grey30"), threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,chr.den.col=c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3"),bin.size=1e6,signal.col=c("red","green"),signal.cex=c(1,1),signal.pch=c(19,19),file="jpg",memo="",dpi=300)

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)