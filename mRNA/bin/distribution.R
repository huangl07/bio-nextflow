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
df<-read.table(opt$input,skip=4,comment.char="=",head = T)
p=ggplot(df)+geom_bar(aes(x=1,y=Tag_count,fill=Group),stat="identity")+theme_classic()+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),)+coord_polar("y")
ggsave(p,file=paste(opt$output,"png",sep="."),device="png")


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
