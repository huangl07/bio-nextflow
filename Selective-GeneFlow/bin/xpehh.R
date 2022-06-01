library(optparse)
options(bitmapType='cairo')
option_list <- list(
  make_option(c("-x", "--xpehh"), type = "character", action = "store", help = "Input xpclr file"),
  make_option(c("-t", "--threshold"), type = "integer", action = "store", default = 0.05, help = "threshold for select")
)
opt = parse_args(OptionParser(option_list = option_list,
                              usage = "Usage: %prog [options] \nDescription: This Script is used to draw xpclr manhattan!"))
if(is.null(opt$xpehh)){opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"),args="--help")}

times<-Sys.time()

library(ggplot2)
library(dplyr)
library(ggpubr)


dfs <- read.table(opt$xpehh, sep="\t",header=TRUE)
colnames(dfs)[1]="chrom"
dfs=dfs %>% filter(xpehh!= "xpehh") %>% filter(xpehh!="")
df=dfs %>% filter(xpehh >quantile(as.numeric(dfs$xpehh),1-opt$threshold,na.rm=T) |  xpehh <quantile(as.numeric(dfs$xpehh),1-opt$threshold,na.rm=T))
pop=gsub(".xpehh.result","",opt$xpehh)
print(pop)
write.csv(file=paste(pop,"xpehh.select",sep="."),df,row.names=F,quote=F)
df=dfs
df$xpehh=as.numeric(df$xpehh)
df$BP=as.numeric(df$pos)/2
quantile1=quantile(df$xpehh,1-opt$threshold,na.rm=T)
quantile2=quantile(df$xpehh,opt$threshold,na.rm=T)

df$CHR=df$chrom
lev=NULL
lev$CHR<-levels(as.factor(df$CHR))
lev$order<-gsub("\\D","",lev$CHR)
lev$order=as.numeric(lev$order)
dfpos=merge(df,lev,by="CHR")
dfpos=arrange(dfpos,order,BP)
dpos <- dfpos %>% group_by(order) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>%
  left_join(dfpos, ., by=c("order"="order")) %>%
  arrange(order, BP) %>%
  mutate( BPcum=BP+tot)
axisdf <- dpos %>% group_by(CHR) %>% summarize(center=(as.numeric(max(BPcum)) + as.numeric(min(BPcum))) / 2 )

  p1 = ggplot(dpos) +geom_point(aes(x=BPcum, y=xpehh,color=as.factor(order)))+geom_hline1(yintercept =quantile,linetype="dashed",col="red")++geom_hline2(yintercept =quantile2,linetype="dashed",col="red")
  p1=p1 + scale_color_manual(values = rep(c("#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2"), nlevels(as.factor(dpos$CHR))))+
		    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
		    theme_bw() +ylab("Dxy")+
		    theme( 
		      legend.position="none",
		      panel.grid.major.x = element_blank(),
		      panel.grid.minor.x = element_blank(),
              axis.title.x=element_blank()
		    )
  pop=gsub(".xpehh.result","",opt$xpehh)
  ggsave(p1,file=paste(pop,"xpehh.manhattan.png",sep="."),device="png",height=3,width=16)
  ggsave(p1,file=paste(pop,"xpehh.manhattan.pdf",sep="."),device="pdf",height=3,width=16)

escaptime=Sys.time()-times
print("Done!")
print(escaptime)

escaptime=Sys.time()-times
print("Done!")
print(escaptime)
