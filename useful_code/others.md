
- [读入数据](#读入数据)
- [对"行"进行操作](#对"行"进行操作)
- [others](#others)
- [统计表格](#统计表格)
- [统计外显子长度](#统计外显子长度)

# 读入数据
```r
# 读入数据记得把抬头和行名处理
data <- read.csv("sample.csv",header=TURE, row.names = 1)
head(data)
```
# 对"行"进行操作
```r
# 去除前5行
countdata <- data[-(1:5),]

# 得到每一行的行名
row_names <- roe.names(countdata)

# 处理行名，去掉后面".w+"
name_replace <- gsub("\\.\\w+", "",row.names(countdata))
# 句号.和加号w+是特殊的，要添加\\来识别
# 表示句号后面的word都要被替换为 [空]
row.names(countdata) <- name_replace

# 去除某些行
countdata <- countdata[rowSums(countdata) > 0,]
```

# others 
```r
1. rm(list=ls())
# 这个命令将从当前的R环境中删除所有对象（变量、函数等）。
# ls()函数用于列出当前环境中的所有对象，list=ls()参数将ls()的结果作为要删除的对象列表传递给rm()函数，从而清空所有对象。

2. getwd("D:/../..")
3. setwd("D:/../..")
# 设置工作目录
```

# 统计表格
每个样本都有一个统计每个基因含有多少个count的表格；有多个样本，想将他们统计起来。
```r
# 得到文件样本编号
files <- list.files(".", "*.count")
f_lists <- list()
for(i in files){
    prefix = gsub("(_\\w+)?\\.count", "", i, perl=TRUE)
    f_lists[[prefix]] = i
}

# 批量导入文件
id_list <- names(f_lists)
data <- list()
count <- 0
for(i in id_list){
  count <- count + 1
  a <- read.table(f_lists[[i]], sep="\t", col.names = c("gene_id",i))
  data[[count]] <- a
}

# 合并文件
data_merge <- data[[1]]
for(i in seq(2, length(id_list))){
    data_merge <- merge(data_merge, data[[i]],by="gene_id")
}

write.csv(data_merge, "merge.csv", quote = FALSE, row.names = FALSE)
```

# 统计外显子长度
1. ref.1
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicFeatures")
library("GenomicFeatures")
## 导入gff3文件
txdb <- makeTxDbFromGFF("genome.gff3",format="gff3")
## 获取外显子位置
exons_gene <- exonsBy(txdb, by = "gene")
## 去除外显子重叠部分，计算外显子长度
gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
## 查看一下
class(gene_lens)
length(gene_lens)
data <- t(as.data.frame(gene_lens))
# 写入文件
write.table(data, file = "mRatBN7.2_gene_len.tsv", row.names = TRUE, sep="\t", quote = FALSE, col.names = FALSE)
```