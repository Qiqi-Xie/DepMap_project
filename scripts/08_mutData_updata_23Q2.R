################################################
################################################
### 作者：果子
### 更新时间：2023-12-06
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人邮箱：hello_guozi@126.com
### 有个问题, BRAF的V600E数据没有

## 数据预处理
## CCLE 究竟如何界定 突变的？
library(dplyr)
library(tidyr)
library(tibble)

#######################################################
### 22Q4V2版本
rm(list = ls())
mutData_raw <- data.table::fread("data/DepMap_Public_20Q4v2/CCLE_mutations.csv",data.table = F)
test <- mutData_raw[mutData_raw$Hugo_Symbol=="BRAF",]
table(mutData_raw$Variant_Classification,mutData_raw$isDeleterious)
table(mutData_raw$Variant_Classification,mutData_raw$Variant_annotation)
table(mutData_raw$Variant_annotation,mutData_raw$isDeleterious)

index <- c("De_novo_Start_OutOfFrame", "Frame_Shift_Del", 
           "Frame_Shift_Ins", "IGR", "Intron", 
           "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", 
           "Start_Codon_Del", "Start_Codon_Ins",
           "Stop_Codon_Del", "Stop_Codon_Ins")

### 1749,19543
mutData <- mutData_raw %>% 
  mutate(isDeleterious = Variant_Classification %in% index) %>% 
  dplyr::select(DepMap_ID,Hugo_Symbol,isDeleterious) %>% 
  mutate(isDeleterious= as.numeric(isDeleterious)) %>%
  group_by(DepMap_ID,Hugo_Symbol) %>% 
  summarise(isDeleterious=sum(isDeleterious)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Hugo_Symbol",
              values_from = "isDeleterious",
              values_fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("DepMap_ID")

saveRDS(mutData,file = "data/DepMap_Public_20Q4v2/mutData_with_Missense.rds")

### 测试
### 修改前,85,704
### 修改后,118,671
### 原文,133,656

mutData <- readRDS(file = "data/DepMap_Public_20Q4v2/mutData.rds")
#mutData <- readRDS(file = "data/DepMap_Public_20Q4v2/mutData_with_Missense.rds")
geneDependency <- readRDS(file = "data/DepMap_Public_20Q4v2/geneDependency.rds")
cellinfor <- readRDS(file = "data/DepMap_Public_20Q4v2/cellinfor.rds")

### 细胞系取交集，789
coID <- intersect(rownames(geneDependency),rownames(cellinfor)) %>% 
  intersect(rownames(mutData))

### 提取数据 ARID1A的突变信息，以及HMGCR的Dependency 信息
data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,"primary_disease"],
  Mutation = ifelse(mutData[coID,"ARID1A"]==0,"WT","Mut"),
  Dependency = geneDependency[coID,"HMGCR"]
)

## 数据筛选，Mut >=2,WT>=5(跟原文不同)
TissueStatus <- data %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,
              values_fill = 0) %>% 
  filter(Mut >=0,WT>=0) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 

#######################################################
### 22Q2版本
rm(list = ls())
mutData_raw <- data.table::fread("data/DepMap_Public_22Q2/CCLE_mutations.csv",data.table = F)
test <- mutData_raw[mutData_raw$Hugo_Symbol=="BRAF",]
table(mutData_raw$Variant_Classification,mutData_raw$isDeleterious)

index <- c("De_novo_Start_OutOfFrame", "Frame_Shift_Del", 
           "Frame_Shift_Ins", "IGR", "Intron", 
           "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", 
           "Start_Codon_Del", "Start_Codon_Ins",
           "Stop_Codon_Del", "Stop_Codon_Ins")

mutData <- mutData_raw %>% 
  mutate(isDeleterious = Variant_Classification %in% index) %>% 
  dplyr::select(DepMap_ID,Hugo_Symbol,isDeleterious) %>% 
  mutate(isDeleterious= as.numeric(isDeleterious)) %>%
  group_by(DepMap_ID,Hugo_Symbol) %>% 
  summarise(isDeleterious=sum(isDeleterious)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Hugo_Symbol",
              values_from = "isDeleterious",
              values_fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("DepMap_ID")

saveRDS(mutData,file = "data/DepMap_Public_22Q2/mutData_with_Missense.rds")

#mutData <- readRDS(file = "data/DepMap_Public_20Q4v2/mutData.rds")
#mutData <- readRDS(file = "data/DepMap_Public_20Q4v2/mutData_with_Missense.rds")
mutData <- readRDS(file = "data/DepMap_Public_22Q2/mutData_with_Missense.rds")
### 22Q4V2
# geneDependency <- readRDS(file = "data/DepMap_Public_20Q4v2/geneDependency.rds")
# cellinfor <- readRDS(file = "data/DepMap_Public_20Q4v2/cellinfor.rds")
### 22Q2
geneDependency <- readRDS(file = "data/DepMap_Public_22Q2/geneDependency.rds")
cellinfor <- readRDS(file = "data/DepMap_Public_22Q2/cellinfor.rds")

### 细胞系取交集，789
coID <- intersect(rownames(geneDependency),rownames(cellinfor)) %>% 
  intersect(rownames(mutData))

### 提取数据 ARID1A的突变信息，以及HMGCR的Dependency 信息
data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,"primary_disease"],
  Mutation = ifelse(mutData[coID,"ARID1A"]==0,"WT","Mut"),
  Dependency = geneDependency[coID,"HMGCR"]
)

## 数据筛选，Mut >=2,WT>=5(跟原文不同)
TissueStatus <- data %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,
              values_fill = 0) %>% 
  filter(Mut >=0,WT>=0) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 

##################################################################
## 新的版本 23Q2 
rm(list = ls())
library(dplyr)
library(tidyr)
library(tibble)
mutData_raw <- data.table::fread("data/DepMap_Public_23Q2/OmicsSomaticMutations.csv",data.table = F)

table(mutData_raw$VariantInfo,mutData_raw$CCLEDeleterious)

test <- mutData_raw[mutData_raw$HugoSymbol=="BRAF",]

index <- c("FRAME_SHIFT_INS", 
           "IN_FRAME_DEL", 
           "MISSENSE", 
           "NONSENSE",
           "NONSTOP", 
           "START_CODON_INS")

## ModelID 就是以前的DepMap_ID
mutData <- mutData_raw %>% 
  mutate(isDeleterious = VariantInfo %in% index) %>% 
  dplyr::select(ModelID,HugoSymbol,isDeleterious) %>% 
  rename(DepMap_ID=ModelID,Hugo_Symbol=HugoSymbol) %>% 
  mutate(isDeleterious= as.numeric(isDeleterious)) %>%
  group_by(DepMap_ID,Hugo_Symbol) %>% 
  summarise(isDeleterious=sum(isDeleterious)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Hugo_Symbol",
              values_from = "isDeleterious",
              values_fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("DepMap_ID")

saveRDS(mutData,file = "data/DepMap_Public_23Q2/mutData.rds")

### 23Q2 有自己的细胞系信息
cellinfor <- data.table::fread(file = "data/DepMap_Public_23Q2/Model.csv",data.table = F)
colnames(cellinfor)[1] = "DepMap_ID"

cellinfor22Q2 <- readRDS(file = "data/DepMap_Public_22Q2/cellinfor.rds")
cellinfor22Q2 <- cellinfor22Q2[,c("DepMap_ID","lineage","primary_disease")]
cellinfor <- merge(cellinfor,cellinfor22Q2,by = "DepMap_ID")
rownames(cellinfor) <- cellinfor$DepMap_ID
saveRDS(cellinfor,file = "data/DepMap_Public_23Q2/cellinfor.rds")

### 也有自己的geneDependency
## 2.DepMap 基因Dependency 数据
rm(list = ls())
### 1095 vs 1086
geneDependency <- data.table::fread("data/DepMap_Public_23Q2/CRISPRGeneDependency.csv",data.table = F)
test <- geneDependency[1:10,1:10]

## 修改列名
colnames(geneDependency) <- gsub("\\s+\\(\\d+\\)","",colnames(geneDependency))
## 第一列变行名
rownames(geneDependency) <- geneDependency[,1]
geneDependency <- geneDependency[,-1]
## 保存数据
saveRDS(geneDependency,file = "data/DepMap_Public_23Q2/geneDependency.rds")

###########################################################################
rm(list = ls())
mutData <- readRDS(file = "data/DepMap_Public_23Q2/mutData.rds")
cellinfor <- readRDS(file = "data/DepMap_Public_23Q2/cellinfor.rds")
geneDependency <- readRDS(file = "data/DepMap_Public_23Q2/geneDependency.rds")

### 细胞系取交集，1072
coID <- intersect(rownames(geneDependency),rownames(cellinfor)) %>% 
  intersect(rownames(mutData))

### 提取数据 ARID1A的突变信息，以及HMGCR的Dependency 信息
### "lineage","primary_disease","OncotreeCode","OncotreeSubtype","OncotreePrimaryDisease","OncotreeLineage"
Tissue = "lineage"
mutGene = "ARID1A"
targetGene = "HMGCR"
data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,Tissue],
  Mutation = ifelse(mutData[coID,mutGene]==0,"WT","Mut"),
  Dependency = geneDependency[coID,targetGene]
)

## 数据筛选，Mut >=3,WT>=5(跟原文不同)
TissueStatus <- data %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,
              values_fill = 0) %>% 
  filter(Mut >=3,WT>=5) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 

mydata <- data %>% 
  filter(Tissue %in% TissueStatus$Tissue)

## 作图数据整理
Mean_Dep <- mydata %>% 
  filter(Mutation=="Mut") %>% 
  group_by(Tissue) %>% 
  summarise(Sensitivity = mean(Dependency))

Diff_Dep <- mydata %>%
  group_by(Tissue) %>%
  summarize(Difference = mean(Dependency[Mutation == "Mut"]) - mean(Dependency[Mutation == "WT"]))

## 数据合并
PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue") %>% 
  merge(TissueStatus,by="Tissue")


## 作图展示
library(ggplot2)
library(ggrepel)
library(stringr)
library(ggfun)

ggplot(PlotData, aes(x = Sensitivity, y = Difference)) +
  geom_point(aes(size = Sum), shape = 21, color = "black", fill = "green", alpha = 0.6, stroke = 1.5) + 
  geom_text_repel(aes(label = Tissue), vjust = 1.5, hjust = 1.5) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "yellow", alpha = 0.1, colour = "green", linetype = "dashed", linewidth = 1) +
  labs(x = "Sensitivity Score for Cell Lines with Mutations", 
       y = "Mutant vs. Wildtype Sensitivity", 
       title = paste0("Mutant (n=", unique(PlotData$Mut_Number),
                      " cell lines) vs. Wildtype (n=",
                      unique(PlotData$WT_Number), " cell lines)")) +
  theme_bw() + 
  guides(size = guide_legend(title = str_wrap("Number of cells", width = 8))) +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        panel.grid.minor = element_line(linetype = "dashed", linewidth = 0.5),
        panel.background = element_blank(),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(4, 4, 4, 4),
        legend.background = element_roundrect(color = "#808080", linetype = 1)
  )

#############################################################
#############################################################
### 数据又更新了！2024年11月23日
### 等待确认
rm(list = ls())
library(dplyr)
library(tidyr)
library(tibble)
mutData_raw <- data.table::fread("data/DepMap_Public_23Q4/OmicsSomaticMutations.csv",data.table = F)
names(table(mutData_raw$VariantInfo))
names(table(mutData_raw$VariantType))

### 推荐阅读
### CCLE 如何标记突变
### https://forum.depmap.org/t/how-is-the-isdeleterious-column-in-the-ccle-mutations-csv-file-determined/129