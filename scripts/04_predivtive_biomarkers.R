################################################
################################################
### 作者：果子
### 更新时间：2023-10-22
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人邮箱：hello_guozi@126.com

### 基因A的Dependdency受到哪些基因调控
### 或者哪些基因的表达能够预测抑制剂疗效

rm(list = ls())
geneEffect <- readRDS(file = "output/depmap_geneEffect_1005_17285.rds")
genedata <- geneEffect[,"ESR1"]

exprSet <- readRDS(file = "output/ccle_exprSet_1005_17285.rds")

cor.test(genedata,exprSet[,"ESR1"],method = "pearson")


corData <- data.frame()
gene1 = "ESR1"
genedata <- geneEffect[,"ESR1"]
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

### 推荐阅读:
### 预测模型构建,lasso,随机森林等机器学习算法
### https://bookdown.org/gongchangzhaojie/TranslationalBioinformaticsWithR/predictive-biomarkers.html