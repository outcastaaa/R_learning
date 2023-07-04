# his 图
```r
HIPP_totaldiff_length0.5 <- read.table("D:/brain/brain/diff_peak/0.5/mouse/HIPP_totaldiff_length.txt")
summary(HIPP_totaldiff_length0.5$V1)
hist(abs(as.numeric(HIPP_totaldiff_length0.5[,1])),\
breaks=500,xlab = "Fragment length(bp)",\
ylab = "Frequency",main = "HIPP_totaldiffpeak length")
```