#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'infile','i',0,'character',
	'outfile','o',0,'character',
	'group','g',0,'character',
	'help','m',0,'logical',
	'list','l',0,'character'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	q(status=1);
}
if ( !is.null(opt$help)) { print_usage(spec) }
if ( is.null(opt$infile) & is.null(opt$list)) { print_usage(spec)}
if ( is.null(opt$outfile)){ print_usage(spec) }

times<-Sys.time()
#
#ldfile<-read.table("pop.stat.gz",head=TRUE)
#ld<-read.table("pop.stat.gz",head=TRUE,na.strings=c("nan","-nan"),comment.char=":")
#files=ldfile$file
#popid=ldfile$popid
#col<-rainbow(5)
#pop.id<-paste("pop",c(1:5))

library(RColorBrewer)
library(ggplot2)
library(dplyr)
color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
col<-color[sample(1:length(color),length(color))]

#for (i in 1:5){
if(is.null(opt$list)){
	ld=read.table(file=opt$infile,head=TRUE,comment.char=":");
	distance=ld$X.Dist
	distance=ld$X.Dist
	R2=ld$Mean_r.2
	n=length(R2)
	HW.st<-c(C=0.01)
	HW.nonlinear<-nls(R2~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
	tt<-summary(HW.nonlinear)
	new.rho<-tt$parameters[1]
	fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
	newLD<-data.frame(distance,fpoints)
	maxld<-max(newLD$fpoints)
	decay05<-newLD$distance[which.min(abs(newLD$fpoints-maxld/2))]
	decay01<-newLD$distance[which.min(abs(newLD$fpoints-0.1))]
	newLD<-newLD[order(newLD$distance),]
	newLD$distance<-newLD$distance/1000
	p=ggplot(newLD)+geom_line(aes(x=distance,y=fpoints))+xlab("Distance(kb)")+ylab("r^2")+ggtitle(label="LD decay",subtitle=paste("decay05",decay05,sep=":"))+theme_bw()+theme(plot.title = element_text(hjust = 0.5))
	ggsave(file=paste(opt$outfile,"png",sep="."),p,width=16,height=9)
	ggsave(file=paste(opt$outfile,"pdf",sep="."),p,width=16,height=9)
}else{
	ldfile=read.table(file=opt$list);
	read<-function(i){
		ld=read.table(file=ldfile[i,2],head=TRUE,comment.char=":");
		distance=ld$X.Dist
		distance=ld$X.Dist
		R2=ld$Mean_r.2
		n=length(R2)
		HW.st<-c(C=0.01)
		HW.nonlinear<-nls(R2~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
		tt<-summary(HW.nonlinear)
		new.rho<-tt$parameters[1]
		fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
		newLD<-data.frame(distance,fpoints)
		newLD<-newLD[order(newLD$distance),]
		newLD$distance<-newLD$distance/1000
		newLD$id=ldfile[i,1]
		newLD
	}
	df=lapply(1:nrow(ldfile),read)
	dfs=do.call(rbind,df)
	ldstat=dfs %>% group_by(id) %>% summarise(decay05=distance[which.min(abs(fpoints-max(fpoints)/2))],decay01=distance[which.min(abs(fpoints-0.1))])
	write.table(ldstat,file="ldstat.xls",sep="\t",quote=F)
	p=ggplot(dfs)+geom_line(aes(x=distance,y=fpoints,col=id))+xlab("Distance(kb)")+ylab("r^2")+ggtitle("LD decay")+theme_bw()+theme(plot.title = element_text(hjust = 0.5))
	ggsave(file=paste(opt$outfile,"png",sep="."),p,width=16,height=9)
	ggsave(file=paste(opt$outfile,"pdf",sep="."),p,width=16,height=9)
}
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
