# usr/bin R

#相关R包载入：
library(ggridges)
library(ggplot2)
library(cols4all)

#本地测试数据读入：
df <- read.csv('iris',header = T)
head(df)