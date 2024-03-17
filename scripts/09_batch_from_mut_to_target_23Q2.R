################################################
################################################
### 作者：果子
### 更新时间：2023-12-06
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人邮箱：hello_guozi@126.com

### 批量找到ARID1A突变后哪个基因受到的影响最大

rm(list = ls())
### 加载数据分析三剑客
library(dplyr)
library(tidyr)
library(tibble)

mutData <- readRDS(file = "data/DepMap_Public_23Q2/mutData.rds")
cellinfor <- readRDS(file = "data/DepMap_Public_23Q2/cellinfor.rds")
geneDependency <- readRDS(file = "data/DepMap_Public_23Q2/geneDependency.rds")
### 细胞系取交集，1072
coID <- intersect(rownames(geneDependency),rownames(cellinfor)) %>% 
  intersect(rownames(mutData))

Tissue = "lineage"
mutGene = "ARID1A"
targetGene = "HMGCR"

# targetGene= "ARID1B"
# targetGene= "SCAP"
# targetGene= "YPEL5"

data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,Tissue],
  Mutation = ifelse(mutData[coID,mutGene]==0,"WT","Mut"),
  Dependency = geneDependency[coID,targetGene]
)

## 数据筛选，Mut >=2,WT>=5(跟原文不同)
TissueStatus <- data %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n(),.groups = "drop") %>% 
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

## 为了批量,定义一个分数
score = mean(PlotData$Sensitivity)*mean(PlotData$Difference)*100

# score = as.numeric(PlotData$Sensitivity %*% PlotData$Sum)*
#   as.numeric(PlotData$Difference %*% PlotData$Sum)/1000
#######################################################
### 写循环
genelist = colnames(geneDependency)
Tissue = "lineage"
mutGene = "ARID1A"

### 增加一个判断
### mutGene 如果 在所有细胞中不突变
### TissueStatus 没有行数，跳过
results = data.frame()
for(i in 1:length(genelist) ){
  print(i)
  targetGene = genelist[i]
  data <- data.frame(
    DepmapID = coID,
    Tissue = cellinfor[coID,Tissue],
    Mutation = ifelse(mutData[coID,mutGene]==0,"WT","Mut"),
    Dependency = geneDependency[coID,targetGene]
  )
  
  ## 数据筛选，Mut >=2,WT>=5(跟原文不同)
  TissueStatus <- data %>% 
    select(Tissue,Mutation) %>%
    group_by(Tissue,Mutation) %>% 
    summarise(n =n(),.groups = "drop") %>% 
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
  
  score = mean(PlotData$Sensitivity)*mean(PlotData$Difference)*100
  # score = as.numeric(PlotData$Sensitivity %*% PlotData$Sum)*
  #   as.numeric(PlotData$Difference %*% PlotData$Sum)/1000
  
  results[i,1] = mutGene
  results[i,2] = targetGene
  results[i,3] = mean(PlotData$Sensitivity)
  results[i,4] = mean(PlotData$Difference)
  results[i,5] = score
}

colnames(results) = c("mutGene","targetGene","mena_Sensitivity","mean_Diff","Score")
saveRDS(results,file = "output/mut2target_results.rds")
#saveRDS(results,file = "output/mut2target_results_dot.rds")

results <- readRDS(file = "output/mut2target_results.rds")

### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3954704/
### GZ07 批量课程