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
library(GAPIT3)
myY <- read.table(opt$trait, head = TRUE)
colnames(myY)[1]="Taxa"
myG=read.table(opt$vcf,head=FALSE)
if(!dir.exists(paste(opt$output,"gBLUP",sep="/"))){dir.create(paste(opt$output,"gBLUP",sep="/"))}
setwd(paste(opt$output,"gBLUP",sep="/"))
myGAPIT1 <- GAPIT(Y=myY[,c(1,2)],G=myG,PCA.total=3,model=c("gBLUP"),file.fragment = 128)
prediction=myGAPIT1$Pred
prediction.ref=prediction[prediction[,3]==1,]
prediction.inf=prediction[prediction[,3]==2,]
YP.ref <- merge(myY, prediction.ref, by.x = "Taxa", by.y = "Taxa")
YP.inf <- merge(myY, prediction.inf, by.x = "Taxa", by.y = "Taxa")
r.ref=cor(as.numeric(as.vector(YP.ref[,2])),as.numeric(as.vector(YP.ref[,6]) ))
r.inf=cor(as.numeric(as.vector(YP.inf[,2])),as.numeric(as.vector(YP.inf[,6]) ))
storage.ref[rep,1]=r.ref
storage.inf[rep,1]=r.inf
storage=cbind(storage.ref,storage.inf)
colnames(storage)=c("Reference","Inference")
write.table(storage, "GAPIT.Cross.Validation.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = NA)
save.RDS(myGAPIT1,file="gBLUP.RDa")

if(!dir.exists(paste(opt$output,"cBLUP",sep="/"))){dir.create(paste(opt$output,"cBLUP",sep="/"))}
setwd(paste(opt$output,"cBLUP",sep="/"))
myGAPIT2 <- GAPIT(Y=myY[,c(1,2)],G=myG,PCA.total=3,model=c("cBLUP"),file.fragment = 128)
prediction=myGAPIT2$Pred
prediction.ref=prediction[prediction[,3]==1,]
prediction.inf=prediction[prediction[,3]==2,]
YP.ref <- merge(myY, prediction.ref, by.x = "Taxa", by.y = "Taxa")
YP.inf <- merge(myY, prediction.inf, by.x = "Taxa", by.y = "Taxa")
r.ref=cor(as.numeric(as.vector(YP.ref[,2])),as.numeric(as.vector(YP.ref[,6]) ))
r.inf=cor(as.numeric(as.vector(YP.inf[,2])),as.numeric(as.vector(YP.inf[,6]) ))
storage.ref[rep,1]=r.ref
storage.inf[rep,1]=r.inf
storage=cbind(storage.ref,storage.inf)
colnames(storage)=c("Reference","Inference")
write.table(storage, "GAPIT.Cross.Validation.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = NA)
save.RDS(myGAPIT2,file="cBLUP.RDa")

if(!dir.exists(paste(opt$output,"sBLUP",sep="/"))){dir.create(paste(opt$output,"sBLUP",sep="/"))}
setwd(paste(opt$output,"sBLUP",sep="/"))
myGAPIT3 <- GAPIT(Y=myY[,c(1,2)],G=myG,PCA.total=3,model=c("sBLUP"),file.fragment = 128)
prediction=myGAPIT2$Pred
prediction.ref=prediction[prediction[,3]==1,]
prediction.inf=prediction[prediction[,3]==2,]
YP.ref <- merge(myY, prediction.ref, by.x = "Taxa", by.y = "Taxa")
YP.inf <- merge(myY, prediction.inf, by.x = "Taxa", by.y = "Taxa")
r.ref=cor(as.numeric(as.vector(YP.ref[,2])),as.numeric(as.vector(YP.ref[,6]) ))
r.inf=cor(as.numeric(as.vector(YP.inf[,2])),as.numeric(as.vector(YP.inf[,6]) ))
storage.ref[rep,1]=r.ref
storage.inf[rep,1]=r.inf
storage=cbind(storage.ref,storage.inf)
colnames(storage)=c("Reference","Inference")
write.table(storage, "GAPIT.Cross.Validation.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = NA)
save.RDS(myGAPIT3,file="sBLUP.RDa")

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)