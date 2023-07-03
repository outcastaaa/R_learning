
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