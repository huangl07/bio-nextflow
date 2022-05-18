library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'map','m',1,'character',
	'loc','l',1,'character',
	'trt','t',1,'character',
	'out','o',1,'character',
	'num','n',1,'character',
	'lod','d',1,'character',
    'btl','b',1,'charater',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript Rqtl_CP.r --map --loc --trt --out --num
	
Usage:
	--map	map file
	--loc	loc file
	--trt	trt file
	--out	out dir
	--num	pm number
    --btl
	--pvalue	select pvalue
	--lod	threshold lod value
	--help		usage
\n")
	q(status=1);
}
times<-Sys.time()
library('qtl');
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$map) ) { print_usage(spec) }
if ( is.null(opt$loc) ) { print_usage(spec) }
if ( is.null(opt$trt) ) { print_usage(spec) }
if ( is.null(opt$num) ) { opt$num=1000; }
if ( is.null(opt$out) ) { opt$out="./";}
if(!dir.exists(opt$out)){dir.create(opt$out)}
if (is.null(opt$method)){opt$method="cim"}
if ( is.null(opt$pvalue) ) { opt$pvalue=0.05}
if(is.null(opt$btl)){opt$btl = F }else{opt$btl=T}
 opt$num=as.numeric(opt$num)
d<-read.cross(mapfile=opt$map,genfile=opt$loc,phefile=opt$trt,format="mapqtl",crosstype="4way")
if(!dir.exists(opt$out)){dir.create(opt$out)}
setwd(opt$out);
d<-jittermap(d)
d<-sim.geno(d)
d<-calc.genoprob(d)
phe.name<-colnames(d$pheno)[2]
pdf(paste(opt$out,"pheno.pdf",sep="."))
plotPheno(d,pheno.col=phe.name)
dev.off()
png(paste(opt$out,"pheno.png",sep="."))
plotPheno(d,pheno.col=phe.name)
dev.off()
model="normal"
if(opt$btl){
    model="binary"
}

scan<-scanone(d,pheno.col=phe.name,model=model);
scan.pm<-scanone(d,pheno.col=phe.name,n.perm=opt$num,model=model);
markerid<-find.marker(d,chr=scan$chr,pos=scan$pos)

pm.result<-summary(scan.pm,alpha=c(0.01,0.05,0.1))
detail=data.frame(marker=markerid,chr=scan$chr,pos=scan$pos,lod=scan$lod,pm1=pm.result[1,1],pm5=pm.result[2,1],pm10=pm.result[3,1])
scan.result<-summary(scan,format="tabByCol",threshold=pm.result[3,1],drop=1)
if(nrow(scan.result$lod) ==0){scan.result<-summary(scan,format="tabByCol",threshold=2,drop=1)}
markerid1<-find.marker(d,scan.result$lod$chr,scan.result$lod$ci.high)
markerid2<-find.marker(d,scan.result$lod$chr,scan.result$lod$ci.low)
qtl<-makeqtl(d,chr=scan.result$lod$chr,pos=scan.result$lod$pos)
fitqtl<-fitqtl(cross=d,qtl=qtl,pheno.col=phe.name,get.est=TRUE)
result=data.frame(scan.result$lod,pos1=markerid1,pos2=markerid2,fitqtl$result.drop)
write.table(detail,file=paste(opt$out,"detail.result",sep="."),row.names=F,quote=F)
write.table(result,file=paste(opt$out,"qtl-result.result",sep="."),row.names=F,quote=F)
write.table(file=paste(phe.name,".pm.csv",sep=""),sep="\t",scan.pm,quote=F,row.names=T);
legend=pm.result
	pdf(file=paste(phe.name,".scan.pdf",sep=""),height=9,width=16)
	plot(scan)
    print("haha")
	abline(h=pm.result,col=rainbow(length(pm.result)))
	legend("topright",legend=legend,col=rainbow(length(pm.result)),pch=1)
	dev.off()
    png(file=paste(phe.name,".scan.png",sep=""),height=9,width=16)
	plot(scan)
	abline(h=pm.result,col=rainbow(length(pm.result)))
	legend("topright",legend=legend,col=rainbow(length(pm.result)),pch=1)
	dev.off()
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
