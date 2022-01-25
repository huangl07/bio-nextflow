times<-Sys.time();
library('getopt');
options(bitmapType='cairo');
spec = matrix(c(
	'matrix','a',0,'character',
	'group','b',0,'character',
   'output','c',0,'character',
	'help' , 'e', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage:
	--matrix	input count matrix file
   --group input group file
   --contrast input condition file [option]
	--output	output dir
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$matrix) ) { print_usage(spec) }
if ( is.null(opt$group) ) { print_usage(spec) }
if ( is.null(opt$output) ) { print_usage(spec) }

library(ggplot2)
library(dplyr)
library(patchwork)
library(pheatmap)
library(Hmisc)
library("RColorBrewer")
group<-read.table(opt$group,head=T,comment.char="^",sep="\t")
gene<-read.csv(opt$matrix,head=T,row.names=1)
if(!dir.exists(opt$output)){dir.create(opt$output)}
setwd(opt$output)
pca<-prcomp(as.matrix(gene))
PCA=as.data.frame(pca$rotation)
PCE=data.frame("PCs"=c(1:length(pca$sdev)),"sdev"=pca$sdev/sum(pca$sdev))
PCA$id=rownames(PCA)
colnames(group)=c("id","group")
PCA=left_join(PCA,group,by="id")
p1=ggplot(PCA)+geom_point(aes(x=PC1,y=PC2,col=group))+stat_ellipse(aes(x=PC1,y=PC2,fill=group),type="norm",geom="polygon",alpha=0.2,color=NA)
p2=ggplot(PCA)+geom_point(aes(x=PC2,y=PC3,col=group))+stat_ellipse(aes(x=PC2,y=PC3,fill=group),type="norm",geom="polygon",alpha=0.2,color=NA)
p3=ggplot(PCA)+geom_point(aes(x=PC1,y=PC3,col=group))+stat_ellipse(aes(x=PC1,y=PC3,fill=group),type="norm",geom="polygon",alpha=0.2,color=NA)
p4=ggplot(PCE)+geom_histogram(aes(x=PCs,y=sdev),stat="identity",fill="lightblue")
p=p1+p2+p3+p4+plot_layout(ncol = 2,nrow=2)
ggsave("PCA.pdf",p,device="pdf",width=16,height=16)
write.csv(file="pca.csv",PCA)
write.csv(file="pce.csv",PCE)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
cor=rcorr(as.matrix(gene))
pdf("cor.pdf")
pheatmap(cor$r,col=colors,cluster_cols=F)
dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()


