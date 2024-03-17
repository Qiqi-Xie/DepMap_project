################################################
################################################
### 作者：果子
### 更新时间：2023-10-22
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人邮箱：hello_guozi@126.com

### Co-Dependency的使用
### top和全部
exprSet <- data.table::fread("data/DepMap_Public_23Q2/CRISPRGeneEffect.csv",data.table = F)
test <- exprSet[1:10,1:10]
cor.test(exprSet[,"SPDEF (25803)"],exprSet[,"ESR1 (2099)"],method = "pearson")

rm(list = ls())
exprSet <- readRDS(file = "output/depmap_geneEffect_1005_17285.rds")
cor.test(exprSet[,"FOXA1"],exprSet[,"ESR1"])
cor.test(exprSet[,"SPDEF"],exprSet[,"ESR1"],method = "pearson")
cor.test(exprSet[,"SPDEF"],exprSet[,"ESR1"],method = "spearman")

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

### 推荐阅读
### 单基因批量相关性分析的妙用
### https://mp.weixin.qq.com/s/TfE2koPhSkFxTWpb7TlGKA
### 单基因批量相关性分析的GSEA
### https://mp.weixin.qq.com/s/sZJPW8OWaLNBiXXrs7UYFw
### 两个功能集成到了GTBAdb中
### http://guotosky.vip:13838/GTBA/

data <- cor(exprSet)

### 拓展使用,新的功能模块:
### A Ubiquitination Cascade Regulating the Integrated Stress Response and Survival in Carcinomas
### A non-canonical tricarboxylic acid cycle underlies cellular identity