
- [读入数据](#读入数据)
- [对"行"进行操作](#对"行"进行操作)
- [others](#others)
- [统计表格](#统计表格)
- [统计外显子长度](#统计外显子长度)
- [PCA_plot](#PCA_plot)
- [PCA_heatmap](#PCA_heatmap)
- [添加表头](#添加表头)
- [组合、计算表格列的内容](组合、计算表格列的内容)

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

# 安装软件包
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager") 
BiocManager::install("DESeq2")
library(DESeq2)
```

# PCA_plot
```r
# 方法1
#构建dds
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ treatment)
#归一化，因为样本量没有超过30，因此用rld
rld <- rlog(dds, blind=FALSE) 
#返回样本名和分组
pcaData <- plotPCA(rld, intgroup=c("treatment"),returnData = T) 
# 按照treatment排序，方便画图标注颜色
pcaData <- pcaData[order(pcaData$treatment,decreasing=F),]
#知道每一个组有多少样本
table(pcaData$treatment)
# control treatment
#      2         6

#根据上面的结果，设置每一个组的数量，方便加颜色；这一步是把每一个样本根据分组情况画出来，效果见图1
# 在横纵坐标写入差异百分比
percentVar <- round(100*attr(pcaData,"percentVar"),1)
plot(pcaData[,1:2],pch = 19,col= c(rep("red",2),rep("green",6)))  
# pcaData[,1:2]代表PCA中的PC1和PC2；
# pch：指定绘制点所使用的符号;
# cex：指定符号的大小。cex是一个数值，表示pch的倍数，默认是1.5倍;
# col：默认绘图颜色。某些函数(如lines、pie)可以接受一个含有颜色值的向量，并自动循环使用。例如：col=c(“red”, “blue”)需要绘制三条线，那么三条颜色分别red、blue、red;
# ylim  y坐标轴的范围，只写出最小值和最大值

#加上样本名字，效果见图2
text(pcaData[,1],pcaData[,2],row.names(pcaData),cex=1, font = 1) 

#加上图例
legend(-70,43,inset = .02,pt.cex= 1.5,title = "Grade",legend = c("treatment", "control"), 
       col = c( "red","green"),pch = 19, cex=0.75,bty="n")



# 方法2
vsdata <- rlog(dds, blind=FALSE)
plotPCA(vsdata, intgroup="condition") + ylim(-20, 20)
text(pcaData[,1],pcaData[,2],row.names(pcaData),cex=0.5, font = 1)
```

# PCA_heatmap
```r
# 对象为dds，分类依据是treatment
# 颜色管理包（不是必须）
library("RColorBrewer")
# 得到数据对象中基因的计数的转化值
gene_data_transform <- assay(rld)
# 使用t()进行转置
# 使用dist方法求样本之间的距离
sampleDists <- dist(t(gene_data_transform))
# 转化为矩阵用于后续pheatmap()方法的输入
sampleDistMatrix <- as.matrix(sampleDists)
# 将矩阵的名称进行修改
rownames(sampleDistMatrix) <- paste(rld$treatment, rld$condition, rld$ids,sep="-")
colnames(sampleDistMatrix) <- paste(rld$treatment, rld$condition, rld$ids,sep="-")
# 没办法显示ids，因为rld没办法导进入，还得想办法搞colnames

# 对象为dds1，分类依据是condition

# 设置色盘
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# 绘制热图与聚类
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

# 合并具有相同列内容的两文件
```r
diff_gene_symbols <- merge(diff_gene_dataframe, rat_symbols, by = c("ensembl_gene_id"))
```

# 写入文件
```r
dir.create("../stat")
write.table(a, "a.tsv", sep="\t", quote = FALSE)
write.table(b, "b.tsv", row.names = F,sep="\t", quote = FALSE)
write.csv(c, "c.csv", row.names = F, quote = FALSE)
```

# MA图
```r
plotMA(result, ylim=c(-10,10))
```

# 添加表头
1. bash
```bash
(echo -e "chr           start            end     methyl%  methyled   unmethyled " && cat ./a.txt)\
 > tem && mv tem a.txt
```
2. r
```r
colnames(rawdata_df) <- c("chr", "start", "end", "methyl%", "methyled", "unmethyled")
```


# 组合、计算表格列的内容
```r
DSS_first_input_data <- first_raw_data %>%
  mutate(chr = paste("chr", chr, sep = "")) %>%
  mutate(pos = start, N = methyled + unmethyled, X = methyled) %>%
  select(chr, pos, N, X)
# %>% 是管道操作符，它将前一个操作的结果作为参数传递给后一个操作
# mutate 函数用于添加、修改或删除数据集的变量（列）
# select(chr, pos, N, X) 选择 chr、pos、N 和 X 这四个变量，生成新的数据集 
```