#!/usr/bin/env Rscript

times<-Sys.time()

library('getopt')
spec<-matrix(c(
    'help','h',0,"logical",
    'best','b',1,"character",
    'best1','1',1,"character",
    'best2','2',1,"character",
    'bestk','k',1,"character",
    'out','o',1,"character"
),byrow=TRUE,ncol=4);
opt=getopt(spec)
print_usage<-function(spec=NULL){
    cat(getopt(spec,usage=TRUE));
    cat("Usage example :\n")
    cat("
Usage example:
    Rscript draw.admixture.r --best --best1 --best2
Usage:
    --best     the best Qvalue
    --best1     the best1 Qvalue     
    --best2     the best2 Qvalue
    --out       output file name
    --bestk     best K file
    --help        usage
\n")
    q(status=1);
}
if(!is.null(opt$help)){print_usage(spec)}
if(is.null(opt$best)){print_usage(spec)}
if(is.null(opt$best1)){opt$best1="best1"}
if(is.null(opt$best2)){opt$best2="best2"}
if(!is.null(opt$bestk)){opt$bestk=as.numeric(opt$bestk)}
library(ggplot2)
library(reshape2)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
Set1 <- brewer.pal(n = 9,name = "Set1")
Set3 <- brewer.pal(n = 12,name = 'Set3')
Set2 <- brewer.pal(n = 8,name = "Set2")
Dark2 <- brewer.pal(n = 8,name = "Dark2")
Paired <- brewer.pal(n = 12,name = "Paired")
Set <- unique(c(Set1,Set3,Set2,Dark2,Paired))
yellow <- c("#FFFF99","#FFED6F","#FFFFB3","#FFFF33")
Set <- Set[!Set %in% yellow]
######Q矩阵文件########
data<-read.table(opt$best,row.names=1)
orders<- do.call(order, as.list(data))
data$id=rownames(data)
df=melt(data)
df$id=factor(df$id,levels=data$id[rev(orders)])
colnames(df)=c("sample","type","Qvalue")
n=1
if(!is.null(opt$best1) & file.exists(opt$best1)){
    data1<-read.table(opt$best1,row.names=1)
    data1$id=rownames(data1)
    df1=melt(data1)
    df1$id=factor(df1$id,levels=levels(df$sample))
    colnames(df1)=c("sample","type","Qvalue")
    n=n+1
}
if(!is.null(opt$best2) & file.exists(opt$best2)){
    data2<-read.table(opt$best2,row.names=1)
    data2$id=rownames(data2)
    df2=melt(data2)
    df2$id=factor(df2$id,levels=levels(df$sample))
    colnames(df2)=c("sample","type","Qvalue")
    n=n+1
}
if(n==1){
    b<-ggplot(df,aes(x=sample,y=Qvalue,fill=type))+
    geom_bar(stat='identity',colour="grey80")+
    scale_fill_manual(values=Set)+
    scale_y_continuous(expand=c(0,0))+ggtitle("Population Structure")+
    theme(legend.position="none",
       plot.title = element_text(hjust = 0.5),
       axis.text.x=element_text(size = 3,angle = 90, vjust = 0.5,
                                   face = "bold", colour = "black", hjust=0.95),
       axis.text.y=element_text(face = "bold", size = 6, colour = "black"),
       axis.title.x=element_blank()
       )
    if(!is.null(opt$bestk)){
        b=b+ylab(paste("k",opt$bestk,sep=" = "))
    }

    ggsave(b,file="admixture.pdf",width=19,height=3)
    ggsave(b,file="admixture.png",width=19,height=3)
}
if(n==2){
    b<-ggplot(df1,aes(x=sample,y=Qvalue,fill=type))+
    geom_bar(stat='identity',colour="grey80")+
    scale_fill_manual(values=Set)+
    scale_y_continuous(expand=c(0,0))+ggtitle("Population Structure")+
    theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
       axis.text.x=element_blank(),
       axis.text.y=element_text( size = 6, colour = "black"),
       axis.title.x=element_blank()
        )
    b1<-ggplot(df,aes(x=sample,y=Qvalue,fill=type))+
    geom_bar(stat='identity',colour="grey80")+
    scale_fill_manual(values=Set)+
    scale_y_continuous(expand=c(0,0))+
    theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
       axis.text.x=element_text(size = 3,angle = 90, vjust = 0.5,
                                   , colour = "black", hjust=0.95),
       axis.text.y=element_text( size = 6, colour = "black"),
       axis.title.x=element_blank()
        )
    if(!is.null(opt$bestk)){
        b=b+ylab(paste("k",opt$bestk-1,sep=" = "));
        b1=b1+ylab(paste("k",opt$bestk,sep=" = "))
    }

    
    p<-b+b1+plot_layout(ncol=1)
    ggsave(p,file="admixture.pdf",width=19,height=6)
    ggsave(p,file="admixture.png",width=19,height=6)
}
if(n==3){
    b<-ggplot(df1,aes(x=sample,y=Qvalue,fill=type))+
    geom_bar(stat='identity',colour="grey80")+
    scale_fill_manual(values=Set)+
    scale_y_continuous(expand=c(0,0))+ggtitle("Population Structure")+
    theme(legend.position="none",
       plot.title = element_text(hjust = 0.5),
       axis.ticks.x= element_blank(),
       axis.text.x=element_blank(),
       axis.text.y=element_text( size = 6, colour = "black"),
       axis.title.x=element_blank()
       )
    b1<-ggplot(df,aes(x=sample,y=Qvalue,fill=type))+
    geom_bar(stat='identity',colour="grey80")+
    scale_fill_manual(values=Set)+
    scale_y_continuous(expand=c(0,0))+
    theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x= element_blank(),
       axis.text.x=element_blank(),
       axis.text.y=element_text( size = 6, colour = "black"),
       axis.title.x=element_blank()
        )
    b2<-ggplot(df2,aes(x=sample,y=Qvalue,fill=type))+
    geom_bar(stat='identity',colour="grey80")+
    scale_fill_manual(values=Set)+
    scale_y_continuous(expand=c(0,0))+
    theme(legend.position="none",
       plot.title = element_text(hjust = 0.5),
       axis.ticks.x= element_blank(),
       axis.text.x=element_text(size = 3,angle = 90, vjust = 0.5,
                                   colour = "black", hjust=0.95),
       axis.text.y=element_text( size = 6, colour = "black"),
       axis.title.x=element_blank()
        )
    if(!is.null(opt$bestk)){
        b=b+ylab(paste("k",opt$bestk-1,sep=" = "));
        b1=b1+ylab(paste("k",opt$bestk,sep=" = "));
        b2=b2+ylab(paste("k",opt$bestk+1,sep=" = "))
    }
    p<-b+b1+b2+plot_layout(ncol=1)
    ggsave(p,file="admixture.pdf",width=19,height=9)
    ggsave(p,file="admixture.png",width=19,height=9)
}

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)