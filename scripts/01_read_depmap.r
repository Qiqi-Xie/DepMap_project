################################################
################################################
### 作者：果子
### 更新时间：2023-10-22
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人邮箱：hello_guozi@126.com

rm(list = ls())
### 把Depmap数据读取进来

## Gene effect
geneEffect <- data.table::fread("data/DepMap_Public_22Q2/CRISPR_gene_effect.csv",data.table = F)
test <- geneEffect[1:10,1:10]

## Gene Dependecny
geneDependency <- data.table::fread("data/DepMap_Public_22Q2/CRISPR_gene_dependency.csv",data.table = F)
test <- geneDependency[1:10,1:10]

## CCLE 细胞系表达量
exprSet <- data.table::fread("data/DepMap_Public_22Q2/CCLE_expression.csv",data.table = F)
test <- exprSet[1:10,1:10]

## CCLE 细胞信息
cellinfor <- data.table::fread("data/DepMap_Public_22Q2/sample_info.csv",data.table = F)

## mutate information
mutData <- data.table::fread("data/DepMap_Public_22Q2/CCLE_mutations.csv",data.table = F)

################################
## 修剪细胞系信息

commonindex <- intersect(geneEffect$DepMap_ID,exprSet$V1)

## geneEffect
rownames(geneEffect) <- geneEffect[,1]
geneEffect <- geneEffect[,-1]
test <- geneEffect[1:10,1:10]
geneEffect <- geneEffect[commonindex,]
colnames(geneEffect) <- gsub("\\s+\\(\\d+\\)","",colnames(geneEffect))

## exprSet
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
exprSet <- exprSet[commonindex,]
colnames(exprSet) <- gsub("\\s+\\(\\d+\\)","",colnames(exprSet))

## cellInfor
rownames(cellinfor) <- cellinfor[,1]
cellinfor <- cellinfor[,-1]
cellinfor <- cellinfor[commonindex,]

## commonGenes
commonGenes <- intersect(colnames(geneEffect),colnames(exprSet))
geneEffect <- geneEffect[,commonGenes]
exprSet <- exprSet[,commonGenes]

## saveRds

saveRDS(geneEffect,file = "output/depmap_geneEffect_1005_17285.rds")
saveRDS(exprSet,file = "output/ccle_exprSet_1005_17285.rds")
saveRDS(cellinfor,file = "output/cellInfor_1005.rds")

