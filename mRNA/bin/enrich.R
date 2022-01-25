times<-Sys.time();
library('getopt');
options(bitmapType='cairo');
spec = matrix(c(
	'list','l',0,'character',
    'anno','g',0,'character',
	'base','b',0,'character',
    'output','o',0,'character',
	'help' , 'e', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage:
    -list <file> different gene
    -anno <file> annotation 
    -base  <database> database
    -output <output> output dir for enrich
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$list) ) { print_usage(spec) }
if ( is.null(opt$anno) ) { print_usage(spec) }
if ( is.null(opt$base) ) {
    opt$base="/mnt/ilustre/users/dna/database/GO/GO.database";
}
if ( is.null(opt$output) ) { print_usage(spec) }

if(!dir.exists(opt$output)){dir.create(opt$output);setwd(opt$output)}
godb<-read.table(opt$base,sep="\t");
annotate<-read.table(opt$anno,sep="\t");
list<-read.table(opt$list,sep="\t");

enrich=enricher(list,TERM2GENE=annotate[,c(2:1)],TERM2NAME=godb);

write.table(enrich,file="enrich.result",sep="\t")
pdf("barplot.pdf")
barplot(enrich, showCategory=30)
dev.off()
pdf("dotplot.pdf")
dotplot(enrich, showCategory=30)
dev.off()
pdf("dotplot.pdf")
dotplot(enrich, showCategory=30)
dev.off()


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()




