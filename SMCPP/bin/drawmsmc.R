#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
options(scipen=200)
spec = matrix(c(
	'infile','i',0,'character',
	'outdir','o',0,'character',
	'year','y',0,'character',
	'help','m',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	q(status=1);
}
if ( !is.null(opt$help)) { print_usage(spec) }
if ( is.null(opt$infile)) { print_usage(spec)}
if ( is.null(opt$outdir)){ print_usage(spec) }
library(tidyvorse)
mu <- 1.25e-8
gen <- 30
list<-read.table(opt$infile);

merge<-function(i){
    df<-read.table(list[i,2],head=T)
    df$pop=list[i,1]
    df
}
lapply(1:nrow(list),merge)

df=lapply(1:nrow(list),merge)
dfs=do.call(rbind,df)
dfs <- dfs %>% 
  mutate(generations = left_time_boundary/mu) %>% 
  mutate(time = generations*gen) %>% 
  mutate(Ne = (1/lambda_00)/mu) 

p=ggplot(dfs,aes(x=generations,y=Ne)) +
  geom_step(data = dfs,size=0.5,aes(group=pop,color=pop),alpha=0.5,linetype=1) + 
  geom_smooth(aes(color=pop,fill=pop),span=0.2,se=FALSE) +  
  scale_x_log10() + 
  theme_pubclean() +
  coord_cartesian(ylim = c(0, 1.4e5),xlim=c(0.8e3,1.2e6), expand = FALSE, clip = "off") +  
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(legend.position = "None", text = element_text(size=12)) + 
  ylab(expression(paste("Effective Population Size, ",N[e]))) +



escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
