#!/usr/bin/Rscript

# 从CPAN下载、加载软件包/库
install.packages("")  #需要引号
install.packages(c("",""))
library() #不需要引号
update.packages()


# 加载数据
# 加载以逗号分割的csv数据
data <- read.csv("data.csv", header = TRUE, row.names = 1, sep = ",")
?read.table()
# 最开始以字符串形式导入数据
data <- read.csv("data.csv", stringasfactors = FALSE)
# 某些列是分类数据，可以变为factor
data$sex <- factor(data$sex)


# 列名、行名赋值
names(data) <- c("sample1","sample2","sample3")


