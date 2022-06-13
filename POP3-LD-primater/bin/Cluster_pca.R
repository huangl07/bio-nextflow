#!/usr/bin/env Rscript
times <- Sys.time()
if (!require("pacman")){
  install.packages("pacman")
}
library(pacman)
pacman::p_load(getopt)
###传参信息
spec <- matrix(c(
  'infile', 'i', 0, 'character',
  'gfile', 'g', 0, 'character',
  'outfile', 'o', 0, 'character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--infile    输入的表格
	--gfile     分组信息表 
	--outfile   输出文件夹
	--help      usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$infile))   { print_usage(spec)}
if ( is.null(opt$gfile))  { print_usage(spec) }
if ( is.null(opt$outfile))  { print_usage(spec) }

pacman::p_load(tidyverse, dendextend, circlize)
if (!dir.exists(opt$outfile)){
  dir.create(opt$outfile)
}

###矩阵对角最大化实现，匈牙利算法，行列数目一样的情况
match_table <- function(table){
  row <- dim(table)[1] #获取行数
  col <- dim(table)[2] #获取列数
  for (i in 1:(row-1)){
    new_table <- table[i:row, i:col]
    maxrow <- rownames(which(new_table == max(new_table), arr.ind = TRUE))[1] #获取最大值的行位置
    maxcol <- colnames(new_table)[which(new_table == max(new_table), arr.ind = TRUE)[1,2]] #获取最大值的列位置
    ##换行
    tmp <- table[i,]
    table[i,] <- table[maxrow, ]
    names(table[i,]) <- maxrow
    table[maxrow, ] <- tmp
    ##换列
    tmp <- table[,i]
    table[,i] <- table[, maxcol]
    names(table[,i]) <- maxcol
    table[, maxcol] <- tmp
  }
  return(table)
}

###导入结果矩阵
pca_result <- read.delim(opt$infile, header = FALSE, row.names = 1)

###标准分组信息表格输入
group_info <- read.csv(opt$gfile, header = FALSE, sep = "\t")
names(group_info) <- c("Sample", "Group")
cluster_n <- length(table(group_info$Group)) #分类类别

##结果表加上分组信息
df1 <- pca_result %>%
  data.frame(., sample_name = rownames(.)) %>% #把行名变成一列
  as_tibble() %>% 
  left_join(group_info, c("sample_name" = "Sample")) %>% 
  select(sample_name, Group)
df1$Group <- as.factor(df1$Group)

###层次聚类
df_d <- dist(pca_result)
idx <- 1:dim(df1)[1]
hc <- hclust(df_d,method = "ward.D")
groups <- str_c("cluster", cutree(hc,k=cluster_n))
final_data <- cbind(df1, cluster=groups)
table <- table(final_data$Group, final_data$cluster)
trans_table <- match_table(table)

###计算分类的准确率
accuracy <- paste0(sum(diag(prop.table(trans_table))) * 100, "%") #判对率

###结果输出
##输出分类准确率
write.table(accuracy, paste0(opt$outfile, "/accuracy.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
##输出聚类结果表
write.table(final_data, paste0(opt$outfile, "/cluster.xls"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

escape_time <- Sys.time()-times
print("Done!")
print(escape_time)
