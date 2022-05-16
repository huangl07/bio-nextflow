#!/usr/bin/env Rscript
# load library
times<-Sys.time()
library(ape)
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
	'outgroup','u',1,"character",
	'raxml', 'r' , 1 , "character"
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
	--raxml	character	the raxml software for draw
	--outgroup	character	the outgroup sample split by ,
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
if ( is.null(opt$group) ){print("none colour tree!");}
library(ggtree)
library(ggplot2)

if (is.null(opt$raxml)){
	raxml<-read.tree(file=opt$infile)
}else{
	raxml<-read.raxml(file=opt$infile)
}
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

if (!is.null(opt$group)){
	gro=read.table(opt$group,header=FALSE);
	col=unique(gro$V2);
	cls<-NULL
	for (i in 1:length(col)){
		cls<-c(cls,list(as.character(gro$V1[gro$V2==col[i]])));
	}
	names(cls)=col
	print(col)
	raxml <- groupOTU(raxml,cls)
	length(cls)
	cls
}
g<-NULL; 
if (!is.null(opt$outgroup)){
	if (!is.null(opt$raxml)){
		raxml@phylo <- root(raxml@phylo, outgroup=which(raxml@phylo$tip.label %in% strsplit(opt$outgroup, split=",")[[1]]))
	}else{
		raxml <- root(raxml, outgroup=which(raxml$tip.label %in% strsplit(opt$outgroup, split=",")[[1]]))
	}
}
g<-ggtree(raxml,size=.2,layout="rectangular")+geom_tiplab(size=.2,linesize=.2,align=TRUE)
if (!is.null(opt$raxml)){
g<-g+geom_text2(size=.5,color="black",aes(subset=!isTip, label=bootstrap,color="black"))
}
if(!is.null(opt$group) && !is.null(opt$outgroup)){
g<-g+aes(color=group)+scale_color_manual(values=coldata[1:(length(cls)+1)])
}
if(!is.null(opt$group) && is.null(opt$outgroup)){
g<-g+aes(color=group)+scale_color_manual(values=coldata[1:(length(cls)+1)])
}


pdf(paste(opt$outfile,".rectangular.tree.pdf",sep=""))
print(g)
dev.off()

png(paste(opt$outfile,".rectangular.tree.png",sep=""))
print(g)
dev.off()
g<-NULL;

g<-ggtree(raxml,size=.2,layout="circular")+geom_tiplab2(size=.2,linesize=.2,aes(angle=angle),align=TRUE)
if (!is.null(opt$raxml)){
g<-g+geom_text2(size=.5,color="black",aes(subset=!isTip, label=bootstrap,color="black"))
}
if(!is.null(opt$group) && !is.null(opt$outgroup)){
g<-g+aes(color=group)+scale_color_manual(values=coldata[1:(length(cls)+1)])
}
if(!is.null(opt$group) && is.null(opt$outgroup)){
g<-g+aes(color=group)+scale_color_manual(values=coldata[1:(length(cls)+1)])
}


pdf(paste(opt$outfile,".circular.tree.pdf",sep=""))
print(g)

dev.off()
png(paste(opt$outfile,".circular.tree.png",sep=""))
print(g)

dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
