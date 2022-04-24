times<-Sys.time();
library('getopt');
options(bitmapType='cairo');
spec = matrix(c(
	'binfile','a',0,'character',
	'output','b',0,'character',
    'popt','p',0,'character',
	'lg','i',0,'character',
	'help' , 'e', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript QC.r --base  --qual  --key  --od
	
Usage:
	--binfile	base distribution file
	--output	base quality file
    --popt      popt name
	--lg 		linkage lg id
	--help		usage
\n")
	q(status=1);
}

if (!is.null(opt$help) ) { print_usage(spec) }
if(is.null(opt$binfile)){print_usage(spec)}
if(is.null(opt$output)){print_usage(spec)}

        library(ASMap)
        library(dplyr)

        df=read.table(opt$binfile,head=T,sep="\t",row.names=1)
        mstdf=mstmap.data.frame(df,pop.type=opt$popt,dist.fun="kosambi",p.value=1)
		rmstmap=mstmap(mstdf)
		print(summary.map(rmstmap))
		mstdf=est.rf(mstdf)
        result=est.map(mstdf)
        rdf=data.frame("id"=gsub("L1.","",names(unlist(result))),"pos"=unlist(result))
		colnames(rdf)=c("group",opt$lg);
        write.table(file=opt$output,rdf,row.names=FALSE,quote=FALSE,sep="\t")
		mstdf=replace.map(mstdf,result)
		names(mstdf$geno)=opt$lg
		write.cross(mstdf,format=c("csvr"),file=opt$lg)

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
