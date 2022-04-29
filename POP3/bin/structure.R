times<-Sys.time()

library('getopt');
options(bitmapType='cairo')

spec = matrix(c(
	'help' , 'h', 0, "logical",
	'best','b',1,"character",
	'best1','1',1,"character",
	'best2','2',1,"character",
	'bestk','k',1,"character",
	'outfile' , 'o', 1, "character"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
Options: 
	--help		NULL 		get this help
	--best 	character 	structure file
	--best1
	--best2
	--outfile 	character 	the filename for output graph [forced]
	\n")
	q(status=1);
}
if(is.null(opt$outfile)){print_usage(spec)}

bestk=as.numeric(opt$bestk)

coldata<-c("#ce7038",
"#6086e7",
"#95b03c",
"#5c58be",
"#d1972c",
"#552f7f",
"#5aa554",
"#bd4ca4",
"#45c097",
"#b075d7",
"#ba9a40",
"#5963ad",
"#9c964b",
"#c881ca",
"#a05f32",
"#66a1e5",
"#bc4538",
"#89295f",
"#ba4758",
"#da6295");#pie(rep(1,n), col=coldata)
tbl<-read.table(opt$best,header=FALSE);
col<-length(colnames(tbl))
data<-tbl[,2:col]
id<-tbl$V1;

library(igraph)
n=1
if(!is.null(opt$best1) & file.exists(opt$best1)){
	tbl1<-read.table(opt$best,header=FALSE);
	data1<-tbl1[,2:col]
	n=n+1
}

if(!is.null(opt$best2) & file.exists(opt$best2)){
	tbl2<-read.table(opt$best,header=FALSE);
	data2<-tbl2[,2:col]
	n=n+1
}


if(n == 1){
	pdf(paste(opt$outfile,".pdf",sep=""),width=16, height=9)
	par(mar=c(max(nchar(as.character(id)))*0.7,5,2,0),mfrow=c(n,1),oma=c(2,2,2,2))
	barplot(t(as.matrix(data)),names.arg=id, col=coldata[1:col-1],main="Individual #", ylab=paste("K=",bestk), border=NA,las=2)
	dev.off()

	png(paste(opt$outfile,".png",sep=""),width=1600, height=900)
	par(mar=c(max(nchar(as.character(id)))*0.7,5,2,0),mfrow=c(n,1),oma=c(2,2,2,2))
	barplot(t(as.matrix(data)),names.arg=id, col=coldata[1:col-1],main="Individual #", ylab=paste("K=",bestk), border=NA,las=2)
	dev.off()
}
if(n==2 & !is.null(opt$best1) & file.exists(opt$best1)){
	pdf(paste(opt$outfile,".pdf",sep=""),width=16, height=9)
	par(mar=c(max(nchar(as.character(id)))*0.7,5,2,0),mfrow=c(n,1),oma=c(2,2,2,2))
	barplot(t(as.matrix(data1)),names.arg=id, col=coldata[1:col-1], main=paste("K=",bestk-1), border=NA,las=2,xaxt= "n")
	barplot(t(as.matrix(data)),names.arg=id, col=coldata[1:col-1], main=paste("K=",bestk), border=NA,las=2)
	dev.off()

	png(paste(opt$outfile,".png",sep=""),width=1600, height=900)
	par(mar=c(max(nchar(as.character(id)))*0.7,5,2,0),mfrow=c(n,1),oma=c(2,2,2,2))
	barplot(t(as.matrix(data1)),names.arg=id, col=coldata[1:col-1],main=paste("K=",bestk-1), border=NA,las=2,xaxt= "n")
	barplot(t(as.matrix(data)),names.arg=id, col=coldata[1:col-1], main=paste("K=",bestk), border=NA,las=2)
	dev.off()
}
if(n == 2 & !is.null(opt$best2) & file.exists(opt$best2)){
	pdf(paste(opt$outfile,".pdf",sep=""),width=16, height=9)
	par(mar=c(max(nchar(as.character(id)))*0.7,5,2,0),mfrow=c(n,1),oma=c(2,2,2,2))
	barplot(t(as.matrix(data)),names.arg=id, col=coldata[1:col-1], main=paste("K=",bestk-1), border=NA,las=2,xaxt= "n")
	barplot(t(as.matrix(data2)),names.arg=id, col=coldata[1:col-1], main=paste("K=",bestk), border=NA,las=2)
	dev.off()

	png(paste(opt$outfile,".png",sep=""),width=1600, height=900)
	par(mar=c(max(nchar(as.character(id)))*0.7,5,2,0),mfrow=c(n,1),oma=c(2,2,2,2))
	barplot(t(as.matrix(data)),names.arg=id, col=coldata[1:col-1],main=paste("K=",bestk-1), border=NA,las=2,xaxt= "n")
	barplot(t(as.matrix(data2)),names.arg=id, col=coldata[1:col-1], main=paste("K=",bestk), border=NA,las=2)
	dev.off()
}
if(n==3){
		pdf(paste(opt$outfile,".pdf",sep=""),width=16, height=9)
	par(mar=c(max(nchar(as.character(id)))*0.7,5,2,0),mfrow=c(n,1),oma=c(2,2,2,2))
	barplot(t(as.matrix(data1)),names.arg=id, col=coldata[1:col-1], main=paste("K=",bestk-1), border=NA,las=2,xaxt= "n")
	barplot(t(as.matrix(data)),names.arg=id, col=coldata[1:col-1], main=paste("K=",bestk), border=NA,las=2,xaxt= "n")
	barplot(t(as.matrix(data2)),names.arg=id, col=coldata[1:col-1], main=paste("K=",bestk+1), border=NA,las=2)
	dev.off()

	png(paste(opt$outfile,".png",sep=""),width=1600, height=900,)
	par(mar=c(max(nchar(as.character(id)))*0.7,5,2,0),mfrow=c(n,1),oma=c(2,2,2,2))
	barplot(t(as.matrix(data1)),names.arg=id, col=coldata[1:col-1],main=paste("K=",bestk-1), border=NA,las=2,xaxt= "n")
	barplot(t(as.matrix(data)),names.arg=id, col=coldata[1:col-1], main=paste("K=",bestk), border=NA,las=2,xaxt= "n")
	barplot(t(as.matrix(data2)),names.arg=id, col=coldata[1:col-1], main=paste("K=",bestk+1), border=NA,las=2)
	dev.off()
}

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
