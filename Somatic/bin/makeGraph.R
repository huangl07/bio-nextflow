#!/usr/bin/env Rscript

#!/usr/bin/env Rscript
times<-Sys.time();
library('getopt');
options(bitmapType='cairo');
spec = matrix(c(
	'ploidy','p',0,'character',
	'ratio','r',0,'character',
	'out','o',0,'character',
	'help' , 'e', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	
Usage:
    --ploidy input cnv file
    --ratio input ratio file
    --out input out file
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$ploidy) ) { opt$ploidy=2 }
if ( is.null(opt$ratio) ) { print_usage(spec) }
if ( is.null(opt$out) ) { print_usage(spec) }

dataTable <-read.table(opt$ratio, header=TRUE);

ratio<-data.frame(dataTable)
ploidy <- as.numeric(opt$ploidy)


png(filename = paste(opt$out,".log2.png",sep = ""), width = 1180, height = 1180,
    units = "px", pointsize = 20, bg = "white", res = NA)
plot(1:10)
op <- par(mfrow = c(5,5))

for (i in c(1:22,'X','Y')) {
	tt <- which(ratio$Chromosome==i)
	if (length(tt)>0) {
	 plot(ratio$Start[tt],log2(ratio$Ratio[tt]),xlab = paste ("position, chr",i),ylab = "normalized copy number profile (log2)",pch = ".",col = colors()[88])
	 tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )
	 points(ratio$Start[tt],log2(ratio$Ratio[tt]),pch = ".",col = colors()[136])
	
	
	tt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)
	 points(ratio$Start[tt],log2(ratio$Ratio[tt]),pch = ".",col = colors()[461])
	 tt <- which(ratio$Chromosome==i)
	 
	 #UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:
	 #points(ratio$Start[tt],log2(ratio$CopyNumber[tt]/ploidy), pch = ".", col = colors()[24],cex=4)
	 
	}
	tt <- which(ratio$Chromosome==i)
	
	#UNCOMMENT HERE TO SEE THE EVALUATED MEDIAN LEVEL PER SEGMENT:
	#points(ratio$Start[tt],log2(ratio$MedianRatio[tt]), pch = ".", col = colors()[463],cex=4)
	
}

dev.off()


png(filename = paste(opt$out,".png",sep = ""), width = 1180, height = 1180,
    units = "px", pointsize = 20, bg = "white", res = NA)
plot(1:10)
op <- par(mfrow = c(5,5))

maxLevelToPlot <- 3
for (i in c(1:length(ratio$Ratio))) {
	if (ratio$Ratio[i]>maxLevelToPlot) {
		ratio$Ratio[i]=maxLevelToPlot;
	}
}


for (i in c(1:22,'X','Y')) {
	tt <- which(ratio$Chromosome==i)
	if (length(tt)>0) {
	 plot(ratio$Start[tt],ratio$Ratio[tt]*ploidy,ylim = c(0,maxLevelToPlot*ploidy),xlab = paste ("position, chr",i),ylab = "normalized copy number profile",pch = ".",col = colors()[88])
	 tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )
	 points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[136])
	
	tt <- which(ratio$Chromosome==i  & ratio$Ratio==maxLevelToPlot & ratio$CopyNumber>ploidy)	
	points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[136],cex=4)
	 
	tt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)
	 points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[461])
	 tt <- which(ratio$Chromosome==i)
	 
	 #UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:
	 #points(ratio$Start[tt],ratio$CopyNumber[tt], pch = ".", col = colors()[24],cex=4)
	 
	}
	tt <- which(ratio$Chromosome==i)
	
	#UNCOMMENT HERE TO SEE THE EVALUATED MEDIAN LEVEL PER SEGMENT:
	#points(ratio$Start[tt],ratio$MedianRatio[tt]*ploidy, pch = ".", col = colors()[463],cex=4)
	
}

dev.off()

pdf(paste(opt$out,".log2.pdf",sep = ""))
plot(1:10)
op <- par(mfrow = c(5,5))

for (i in c(1:22,'X','Y')) {
	tt <- which(ratio$Chromosome==i)
	if (length(tt)>0) {
	 plot(ratio$Start[tt],log2(ratio$Ratio[tt]),xlab = paste ("position, chr",i),ylab = "normalized copy number profile (log2)",pch = ".",col = colors()[88])
	 tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )
	 points(ratio$Start[tt],log2(ratio$Ratio[tt]),pch = ".",col = colors()[136])
	
	
	tt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)
	 points(ratio$Start[tt],log2(ratio$Ratio[tt]),pch = ".",col = colors()[461])
	 tt <- which(ratio$Chromosome==i)
	 
	 #UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:
	 #points(ratio$Start[tt],log2(ratio$CopyNumber[tt]/ploidy), pch = ".", col = colors()[24],cex=4)
	 
	}
	tt <- which(ratio$Chromosome==i)
	
	#UNCOMMENT HERE TO SEE THE EVALUATED MEDIAN LEVEL PER SEGMENT:
	#points(ratio$Start[tt],log2(ratio$MedianRatio[tt]), pch = ".", col = colors()[463],cex=4)
	
}

dev.off()


pdf(paste(opt$out,".log2.pdf",sep = ""))
plot(1:10)
op <- par(mfrow = c(5,5))

maxLevelToPlot <- 3
for (i in c(1:length(ratio$Ratio))) {
	if (ratio$Ratio[i]>maxLevelToPlot) {
		ratio$Ratio[i]=maxLevelToPlot;
	}
}


for (i in c(1:22,'X','Y')) {
	tt <- which(ratio$Chromosome==i)
	if (length(tt)>0) {
	 plot(ratio$Start[tt],ratio$Ratio[tt]*ploidy,ylim = c(0,maxLevelToPlot*ploidy),xlab = paste ("position, chr",i),ylab = "normalized copy number profile",pch = ".",col = colors()[88])
	 tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )
	 points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[136])
	
	tt <- which(ratio$Chromosome==i  & ratio$Ratio==maxLevelToPlot & ratio$CopyNumber>ploidy)	
	points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[136],cex=4)
	 
	tt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)
	 points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[461])
	 tt <- which(ratio$Chromosome==i)
	 
	 #UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:
	 #points(ratio$Start[tt],ratio$CopyNumber[tt], pch = ".", col = colors()[24],cex=4)
	 
	}
	tt <- which(ratio$Chromosome==i)
	
	#UNCOMMENT HERE TO SEE THE EVALUATED MEDIAN LEVEL PER SEGMENT:
	#points(ratio$Start[tt],ratio$MedianRatio[tt]*ploidy, pch = ".", col = colors()[463],cex=4)
	
}

dev.off()



