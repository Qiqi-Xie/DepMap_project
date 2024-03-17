################################################
################################################
### 作者：果子
### 更新时间：2023-10-22
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人邮箱：hello_guozi@126.com

### 基因A的表达量调控哪些基因的Dependdency
rm(list = ls())
exprSet <- readRDS(file = "output/ccle_exprSet_1005_17285.rds")
genedata <- exprSet[,"ESR1"]

geneEffect <- readRDS(file = "output/depmap_geneEffect_1005_17285.rds")

corData <- data.frame()
gene1 = "ESR1"
genedata <- exprSet[,"ESR1"]
for (i in 1:ncol(exprSet)) {
  
  ## 1指示
  print(i)
  
  ## 2计算
  gene2 = colnames(geneEffect)[i]
  dd = cor.test(genedata,geneEffect[,gene2])
  
  ## 3.储存
  corData[i,1] = gene1
  corData[i,2] = gene2
  corData[i,3] = dd$estimate
  corData[i,4] = dd$p.value
}

colnames(corData) <- c("Gene1","Gene2","cor","pvalue")

### GSEA
library(clusterProfiler)
mygeneList <- -corData$cor
names(mygeneList) <- corData$Gene2
mygeneList <- sort(mygeneList,decreasing = T)
head(mygeneList)

geneSet <- read.gmt("resource/geneSets/h.all.v2023.2.Hs.symbols.gmt")

mygsea <- GSEA(geneList = mygeneList,TERM2GENE = geneSet)
data <- as.data.frame(mygsea)
library(ggplot2)
dotplot(mygsea,showCategory=30,
        split=".sign",
        font.size = 8,
        label_format = 60)+facet_grid(~.sign)

library(enrichplot)
gseaplot2(mygsea,"HALLMARK_MYC_TARGETS_V2",color = "red",pvalue_table = T)

########################################################
### IL6R
rm(list = ls())
exprSet <- readRDS(file = "output/ccle_exprSet_1005_17285.rds")
genedata <- exprSet[,"IL6R"]

geneEffect <- readRDS(file = "output/depmap_geneEffect_1005_17285.rds")

corData <- data.frame()
gene1 = "IL6R"
genedata <- exprSet[,"IL6R"]
for (i in 1:ncol(exprSet)) {
  
  ## 1指示
  print(i)
  
  ## 2计算
  gene2 = colnames(geneEffect)[i]
  dd = cor.test(genedata,geneEffect[,gene2])
  
  ## 3.储存
  corData[i,1] = gene1
  corData[i,2] = gene2
  corData[i,3] = dd$estimate
  corData[i,4] = dd$p.value
}

colnames(corData) <- c("exp","Dependency","cor","pvalue")

### GSEA
library(clusterProfiler)
mygeneList <- -corData$cor
names(mygeneList) <- corData$Dependency
mygeneList <- sort(mygeneList,decreasing = T)
head(mygeneList)

geneSet <- read.gmt("resource/geneSets/c2.cp.reactome.v2023.2.Hs.symbols.gmt")

mygsea <- GSEA(geneList = mygeneList,TERM2GENE = geneSet)
data <- as.data.frame(mygsea)
library(ggplot2)
dotplot(mygsea,showCategory=30,
        split=".sign",
        font.size = 8,
        label_format = 60)+facet_grid(~.sign)

library(enrichplot)
gseaplot2(mygsea,"REACTOME_CHROMOSOME_MAINTENANCE",color = "red",pvalue_table = T)

### 推荐阅读:
### cGAS–STING drives the IL-6-dependent survival of chromosomally instable cancers
### https://www.bioconductor.org/packages/release/bioc/vignettes/decoupleR/inst/doc/pw_bk.html
### https://saezlab.github.io/progeny/articles/progeny.html

### 拓展方向
### 机器学习
### 协同致死
### 药物敏感性
### 基因功能模块