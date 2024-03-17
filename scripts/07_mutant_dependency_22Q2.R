################################################
################################################
### 作者：果子
### 更新时间：2023-11-22
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人邮箱：hello_guozi@126.com

############################################################
## 数据预处理
## 1.突变信息
rm(list = ls())
mutData_raw <- data.table::fread("data/DepMap_Public_22Q2/CCLE_mutations.csv",data.table = F)
library(dplyr)
library(tidyr)
library(tibble)
mutData <- mutData_raw %>% 
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

saveRDS(mutData,file = "data/DepMap_Public_22Q2/mutData.rds")

## 2.DepMap 基因Dependency 数据
rm(list = ls())
geneDependency <- data.table::fread("data/DepMap_Public_22Q2/CRISPR_gene_dependency.csv",data.table = F)
test <- geneDependency[1:10,1:10]
## 修改列名
colnames(geneDependency) <- gsub("\\s+\\(\\d+\\)","",colnames(geneDependency))
## 第一列变行名
rownames(geneDependency) <- geneDependency[,1]
geneDependency <- geneDependency[,-1]
## 保存数据
saveRDS(geneDependency,file = "data/DepMap_Public_22Q2/geneDependency.rds")


## 3.细胞系信息
rm(list = ls())
cellinfor <- data.table::fread("data/DepMap_Public_22Q2/sample_info.csv",data.table = F)
rownames(cellinfor) <- cellinfor$DepMap_ID

saveRDS(cellinfor,file = "data/DepMap_Public_22Q2/cellinfor.rds")

###############################################################
###############################################################
### 图复现: ARID1A 和HMGCR 的例子
rm(list = ls())

### 加载R包和数据
library(dplyr)
library(tidyr)
library(tibble)
mutData <- readRDS(file = "data/DepMap_Public_22Q2/mutData.rds")
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
  filter(Mut >=2,WT>=5) %>%
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
        panel.grid.minor = element_line(linetype = "dashed", size = 0.5),
        panel.background = element_blank(),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(4, 4, 4, 4),
        legend.background = element_roundrect(color = "#808080", linetype = 1)
  )


##############################################################
### 换个基因 HMGCS1

data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,"primary_disease"],
  Mutation = ifelse(mutData[coID,"ARID1A"]==0,"WT","Mut"),
  Dependency = geneDependency[coID,"HMGCS1"]
)

## 数据筛选，Mut >=2,WT>=5(跟原文不同)
TissueStatus <- data %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,values_fill = 0) %>% 
  filter(Mut >=2,WT>=5) %>%
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
PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue")
PlotData <- merge(PlotData,TissueStatus,by="Tissue")

## 作图展示
library(ggplot2)
library(ggrepel)
library(stringr)

ggplot(PlotData, aes(x = Sensitivity, y = Difference)) +
  geom_point(aes(size = Sum), shape = 21, color = "black", fill = "green", alpha = 0.6, stroke = 1.5) + 
  geom_text_repel(aes(label = Tissue), vjust = 1.5, hjust = 1.5) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "yellow", alpha = 0.1, colour = "green", linetype = "dashed", size = 1) +
  labs(x = "Sensitivity Score for Cell Lines with Mutations", 
       y = "Mutant vs. Wildtype Sensitivity", 
       title = paste0("Mutant (n=", unique(PlotData$Mut_Number),
                      " cell lines) vs. Wildtype (n=",
                      unique(PlotData$WT_Number), " cell lines)")) +
  theme_bw() + 
  guides(size = guide_legend(title = str_wrap("Number of cells", width = 8))) +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        panel.grid.minor = element_line(linetype = "dashed", size = 0.5),
        panel.background = element_blank(),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(4, 4, 4, 4),
        legend.background = element_roundrect(color = "#808080", linetype = 1)
  )


###################################################
### 提取典型组织单独画图

library(ggplot2)
df = mydata[mydata$Tissue=="Endometrial/Uterine Cancer",]
df$Mutation <- factor(df$Mutation, levels = c("WT", "Mut"),
                      labels = c("ARID1A(wt) \n (n=11)", "ARID1A(mut) \n (n=15)"))
ggplot(df, aes(x = Mutation, y = Dependency, fill = Mutation)) +
  geom_violin(alpha = 0.5) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5) +
  labs(y = "Dependency Score")+
  theme_bw() +
  theme(legend.position = "none")

###################################################
### 写成函数

mutPlot <- function(mutGene,targerGene){
  data <- data.frame(
    DepmapID = coID,
    Tissue = cellinfor[coID,"primary_disease"],
    Mutation = ifelse(mutData[coID,mutGene]==0,"WT","Mut"),
    Dependency = geneDependency[coID,targerGene]
  )
  
  ## 数据筛选，Mut >=2,WT>=5(跟原文不同)
  TissueStatus <- data %>% 
    select(Tissue,Mutation) %>%
    group_by(Tissue,Mutation) %>% 
    summarise(n =n()) %>% 
    ungroup() %>% 
    pivot_wider(names_from = "Mutation",
                values_from = n,values_fill = 0) %>% 
    filter(Mut >=2,WT>=5) %>%
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
  PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue")
  PlotData <- merge(PlotData,TissueStatus,by="Tissue")
  
  ## 作图展示
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  
  ggplot(PlotData, aes(x = Sensitivity, y = Difference)) +
    geom_point(aes(size = Sum), shape = 21, color = "black", fill = "green", alpha = 0.6, stroke = 1.5) + 
    geom_text_repel(aes(label = Tissue), vjust = 1.5, hjust = 1.5) + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
             fill = "yellow", alpha = 0.1, colour = "green", linetype = "dashed", size = 1) +
    labs(x = "Sensitivity Score for Cell Lines with Mutations", 
         y = "Mutant vs. Wildtype Sensitivity", 
         title = paste0("Mutant (n=", unique(PlotData$Mut_Number),
                        " cell lines) vs. Wildtype (n=",
                        unique(PlotData$WT_Number), " cell lines)")) +
    theme_bw() + 
    guides(size = guide_legend(title = str_wrap("Number of cells", width = 8))) +
    theme(panel.grid.major = element_line(linetype = "dashed"),
          panel.grid.minor = element_line(linetype = "dashed", size = 0.5),
          panel.background = element_blank(),
          legend.position = c(.95, .05),
          legend.justification = c("right", "bottom"),
          legend.box.just = "right",
          legend.margin = margin(4, 4, 4, 4),
          legend.background = element_roundrect(color = "#808080", linetype = 1)
    )
  
} 


mutPlot(mutGene = "ARID1A",targerGene = "HMGCR")
mutPlot(mutGene = "TP53",targerGene = "FOXA1")

### 能不能批量提取？
### GZ07
### 1. 已知突变基因，找靶基因
### 2. 已知感兴趣的基因，找哪个基因突变后对他影响最大
### 批量的指标是什么？

