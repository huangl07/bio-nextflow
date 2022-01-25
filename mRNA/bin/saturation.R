times<-Sys.time();
library('getopt');
options(bitmapType='cairo');
spec = matrix(c(
	'input','a',0,'character',
	'output','b',0,'character',
	'help' , 'e', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage:
	--input	input epkm.file
	--output	output	draw file
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$input) ) { print_usage(spec) }
if ( is.null(opt$output) ) { print_usage(spec) }
library(ggplot2)
library(dplyr)
library(reshape2)
df<-read.table(opt$input,head=T,comment.char="^")
df=df[,c(4,7:26)]
colnames(df)=c("name",seq(5,100,5))
df$group=as.character(cut(df[,21],breaks=c(0,0.3,0.6,3.5,15,60),include.lowest=F,right=F))
df$group[is.na(df$group)]=">60"
df = df[df[,21] > 0,]
df[,2:21]=df[,2:21]/df[,21]
mdf1=melt(df) %>% group_by(group,variable) %>% filter(abs(value -1) < 0.15) %>% summarise(value=n())
mdf2=df %>% group_by(group) %>% summarise(value=n())
mdf=left_join(mdf1,mdf2,by="group")
mdf$value=mdf$value.x/mdf$value.y
p=ggplot(mdf)+geom_line(aes(x=variable,y=value,group=group,col=group))+coord_polar("y")
ggsave(file=paste(opt$output,"pdf",sep="."),device="pdf",p)
ggsave(file=paste(opt$output,"png",sep="."),device="png",p)

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
