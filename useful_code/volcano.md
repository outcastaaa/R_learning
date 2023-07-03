# 火山图
```r
data <- read.table("volcano.txt",head=F)
pvalue <- -log10(data$V3)
data$chang <- ifelse(abs(data$V2) >0.5 & data$V3 < 0.05, ifelse(data$V2 >0.5, "up","down"), "no")
ggplot(data, aes(x=V2, y=pvalue,color= chang)) +
  geom_point(alpha=.5) +
  theme(panel.grid.major = element_blank(),
        axis.ticks.x = element_line(size=1),
        axis.text.x = element_text(angle=30, hjust=1, vjust=1),
        axis.title.x=element_text(face="italic", colour="darkred", size=14), # 字体
        axis.line = element_line(color="black"),
        plot.title = element_text(colour="red", size=8, face="bold")) +
  ylab("-log10(padj)") + 
  xlab("log2FC") +
  ggtitle("differencial genes") +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = -0.5:0.5) + 
  scale_color_manual(values = c("blue", "grey", "red"))
```