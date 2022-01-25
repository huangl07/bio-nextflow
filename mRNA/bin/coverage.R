times<-Sys.time();
library('getopt');
options(bitmapType='cairo');
spec = matrix(c(
	'file','a',0,'character',
    'output','c',0,'character',
 	'help' , 'e', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage:
	--file	input count matrix file
 	--output	output dir
	--help		usage
\n")
	q(status=1);
}


file<-read.table(opt$file)
read<-function(i){
    df<-read.table(file[i,],head=T)
    df1=data.frame(t(df[,-1]))
    colnames(df1)=df[1,1]
    df1
}
res=sapply(c(1:length(file)),read,simplify=T)
res$percentage=c(1:100)
draw=melt(res,id.vars=c("percentage"))
draw$variable=gsub(".sorted","",draw$variable)
p=ggplot2(draw)+geom_line(aes(x=percentage,y=value,col=variable))+labs(col="sampleID")
ggsave(paste(opt$output,"pdf",sep="."),device="pdf",p)
ggsave(paste(opt$output,"png",sep="."),device="png",p)

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()

