################################################
################################################
### 作者：果子
### 更新时间：2023-10-22
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人邮箱：hello_guozi@126.com

### 复习批量操作
### 基因和基因的相关性

rm(list = ls())

exprSet <- readRDS(file = "output/ccle_exprSet_1005_17285.rds")
cor.test(exprSet[,"FOXA1"],exprSet[,"ESR1"])

dd <- cor.test(exprSet[,"FOXA1"],exprSet[,"ESR1"])
dd$p.value
dd$estimate

gene <- "ESR1"

corData <- data.frame()
gene1 = "ESR1"
genedata <- exprSet[,gene1]

for (i in 1:ncol(exprSet)) {
  
  ## 1指示
  print(i)
  
  ## 2计算
  gene2 = colnames(exprSet)[i]
  dd = cor.test(genedata,exprSet[,gene2])
  
  ## 3.储存
  corData[i,1] = gene1
  corData[i,2] = gene2
  corData[i,3] = dd$estimate
  corData[i,4] = dd$p.value
}

colnames(corData) <- c("Gene1","Gene2","cor","pvalue")

### 推荐阅读: GZ07,批量技能