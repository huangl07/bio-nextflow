#!/usr/bin/env Rscript 
times<-Sys.time()
library("scatterplot3d")
library('getopt');
options(bitmapType='cairo')
#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help' , 'h', 0, "logical",
	'infile' , 'i', 1, "character",
	'outfile' , 'o', 1, "character",
	'group', 'g', 1 , "character",
	'varfile', 'v', 1 , "character"
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
	--infile 	character 	the input file [forced]
	--outfile 	character 	the filename for output graph [forced]
	--group		character	the group file for draw
	--varfile	character	the val file
	\n")
	q(status=1);
}
# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }
# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }
#if ( is.null(opt$group) )	{ print_usage(spec) }
if ( is.null(opt$varfile) )	{ print_usage(spec) }
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
"#da6295");
value<-read.table(opt$varfile,header=FALSE);
value<-value[0:20,]
sum<-sum(value)
value<-value/sum*100
pdf(paste(opt$outfile,".val.pdf",sep=""),width=800,height=800);
valmax=length(value);
barplot(value,col="blue",xlab="PCAs",ylab="Variance%",beside=FALSE,names.arg=c(1:valmax),border=TRUE,ylim=c(0,100),main="Variance of 
PCAs")
dev.off()
value<-signif(value,4)

pca.file<-read.table(opt$infile,header=FALSE)
colnames(pca.file)=c("FID","IID",paste("PC",rep(1:(length(pca.file[1,])-2)),sep=""))
if (!is.null(opt$group)){
	pop.list<-read.table(opt$group,head=FALSE)
	names(pop.list)=c("id","popid")
	pop.id<-as.vector(unique(pop.list$popid));
	color.list<-coldata[1:length(pop.id)];
	for (i in 1:length(pop.id)){
		pop.list$colour[pop.list$popid==pop.id[i]]=color.list[i]
	}
	pop.list=pop.list[order(pop.list[,1]),]
	print(pop.list)
}else{
	pop.id<-c(1);
	color.list<-rainbow(1);
}

pca.file=pca.file[order(pca.file[,1]),]
pdf(paste(opt$outfile,".pc1vspc2.pdf",sep=""))
if (length(pop.id) >1){
	plot(pca.file$PC1,pca.file$PC2,col=pop.list$colour,xlab=paste("PC1(",value[1],"%)",sep=""),ylab=paste("PC2(",value[2],"%)",sep=""),main="PC1 vs PC2")
	legend("right",col=color.list,legend=pop.id,bty=T,pch=1,cex=0.5)
}else{
	plot(pca.file$PC1,pca.file$PC2,xlab=paste("PC1(",value[1],"%)",sep=""),ylab=paste("PC2(",value[2],"%)",sep=""),main="PC1 vs PC2")
}
dev.off();
png(paste(opt$outfile,".pc1vspc2.png",sep=""))
if (length(pop.id) >1){
	plot(pca.file$PC1,pca.file$PC2,col=pop.list$colour,xlab=paste("PC1(",value[1],"%)",sep=""),ylab=paste("PC2(",value[2],"%)"),main="PC1 vs PC2")
	legend("right",col=color.list,legend=pop.id,bty=T,pch=1,cex=0.5)
}else{
	plot(pca.file$PC1,pca.file$PC2,xlab=paste("PC1(",value[1],"%)",sep=""),ylab=paste("PC2(",value[2],"%)",sep=""),main="PC1 vs PC2")
}
dev.off();
pdf(paste(opt$outfile,".pc1vspc3.pdf",sep=""))
if (length(pop.id) >1){
	plot(pca.file$PC1,pca.file$PC3,col=pop.list$colour,xlab=paste("PC1(",value[1],"%)",sep=""),ylab=paste("PC3(",value[3],"%)",sep=""),main="PC1 vs PC2")
	legend("right",col=color.list,legend=pop.id,bty=T,pch=1,cex=0.5)
}else{
	plot(pca.file$PC1,pca.file$PC3,xlab=paste("PC1(",value[1],"%)",sep=""),ylab=paste("PC3(",value[3],"%)",sep=""),main="PC1 vs PC3")
}
dev.off();
png(paste(opt$outfile,".pc1vspc3.png",sep=""))
if (length(pop.id) >1){
	plot(pca.file$PC1,pca.file$PC3,col=pop.list$colour,xlab=paste("PC1(",value[1],"%)",sep=""),ylab=paste("PC3(",value[3],"%)",sep=""),main="PC1 vs PC3")
	legend("right",col=color.list,legend=pop.id,bty=T,pch=1,cex=0.5)
}else{
	plot(pca.file$PC1,pca.file$PC3,xlab=paste("PC1(",value[1],"%)",sep=""),ylab=paste("PC3(",value[3],"%)",sep=""),main="PC1 vs PC3")
}
dev.off();
pdf(paste(opt$outfile,".pc2vspc3.pdf",sep=""))
if (length(pop.id) >1){
	plot(pca.file$PC2,pca.file$PC3,col=pop.list$colour,xlab=paste("PC2(",value[2],"%)",sep=""),ylab=paste("PC3(",value[3],"%)",sep=""),main="PC2 vs PC3")
	legend("right",col=color.list,legend=pop.id,bty=T,pch=1,cex=0.5)
}else{
	plot(pca.file$PC2,pca.file$PC3,xlab=paste("PC2(",value[2],"%)",sep=""),ylab=paste("PC3(",value[3],"%)",sep=""),main="PC2 vs PC3")
}
dev.off()
png(paste(opt$outfile,".pc2vspc3.png",sep=""))
if (length(pop.id) >1){
	plot(pca.file$PC2,pca.file$PC3,col=pop.list$colour,xlab=paste("PC2(",value[2],"%)",sep=""),ylab=paste("PC3(",value[3],"%)",sep=""),main="PC2 vs PC3")
	legend("right",col=color.list,legend=pop.id,bty=T,pch=1,cex=0.5)
}else{
	plot(pca.file$PC2,pca.file$PC3,xlab=paste("PC2(",value[2],"%)",sep=""),ylab=paste("PC3(",value[3],"%)",sep=""),main="PC2 vs PC3")
}
dev.off();
png(paste(opt$outfile,".val.png",sep=""),width=800,height=800);
barplot(value,col="blue",xlab="PCAs",ylab="Variance%",beside=FALSE,names.arg=c(1:valmax),border=TRUE,ylim=c(0,100),main="Variance of 
PCAs")
dev.off()
pdf(paste(opt$outfile,".3D.pdf",sep=""))
if (length(pop.id) >1){
	scatterplot3d(pca.file$PC1,pca.file$PC2,pca.file$PC3,xlab=paste("PC1(",value[1],"%)",sep=""),ylab=paste("PC2(",value[2],"%)",sep=""),zlab=paste("PC3(",value[3],"%)",sep="") ,pch=20,color=pop.list$colour,cex=1,cex.lab=1.4, cex.axis=1.2,lwd=3,angle=55,scale.y=0.7)
	 legend("right",col=color.list,legend=pop.id,bty=T,pch=1,cex=0.5)
}else{																			
	scatterplot3d(pca.file$PC1,pca.file$PC2,pca.file$PC3,xlab=paste("PC1(",value[1],"%)",sep=""),ylab=paste("PC2(",value[2],"%)",sep=""),zlab=paste("PC3(",value[3],"%)",sep="") ,pch=20,cex=1,cex.lab=1.4, cex.axis=1.2,lwd=3,angle=55,scale.y=0.7)
}																			
dev.off()																		
png(paste(opt$outfile,".3D.png",sep=""))														
if (length(pop.id) >1){																	
	scatterplot3d(pca.file$PC1,pca.file$PC2,pca.file$PC3,xlab=paste("PC1(",value[1],"%)",sep=""),ylab=paste("PC2(",value[2],"%)",sep=""),zlab=paste("PC3(",value[3],"%)",sep="") ,pch=20,color=pop.list$colour,cex=1,cex.lab=1.4, cex.axis=1.2,lwd=3,angle=55,scale.y=0.7)
	 legend("right",col=color.list,legend=pop.id,bty=T,pch=1,cex=0.5)
}else{																			
	scatterplot3d(pca.file$PC1,pca.file$PC2,pca.file$PC3,xlab=paste("PC1(",value[1],"%)",sep=""),ylab=paste("PC2(",value[2],"%)",sep=""),zlab=paste("PC3(",value[3],"%)",sep="") ,pch=20,cex=1,cex.lab=1.4, cex.axis=1.2,lwd=3,angle=55,scale.y=0.7)
}
dev.off()


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
