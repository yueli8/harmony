#不建議使用,出來的圖不漂亮,最後出來的結果,不能代入我已有的代碼.
#可以節約integrate的運行時間

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(dplyr)
library(Matrix)
#library(muscat)
library(reshape2)
library(celldex)
library(clustree)
library(SingleR)
setwd("E:/Rstudio_default_working/gse") 

#先構建文件夾,放每一套數據.
library(harmony)
rm(list=ls())
dir=c('TIL.Pt1','TIL.Pt2','TIL.Pt3','TIL.Pt4','TIL.Pt5','TIL.Pt6','TIL.Pt7','TIL.Pt8','PBL.Pt1',
      'PBL.Pt2','PBL.Pt3','PBL.Pt4','PBL.Pt5','PBL.Pt6','PBL.Pt7','PBL.Pt8')
sample_name <- c("blood_P1","blood_P2","blood_P3","blood_P4","blood_P5","blood_P6","blood_P7","blood_P8",
                 "tumor_T1","tumor_T2","tumor_T3","tumor_T4","tumor_T5","tumor_T6","tumor_T7","tumor_T8")
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=sample_name[i], min.cells=3, min.features = 200)
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = percent.mt < 10) }

saveRDS(scRNAlist, "scRNAlist.rds")
library(harmony)
scRNAlist <- readRDS("scRNAlist.rds")
scRNA_harmony <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]], scRNAlist[[4]], scRNAlist[[5]],scRNAlist[[6]], scRNAlist[[7]], 
                                           scRNAlist[[8]], scRNAlist[[9]], scRNAlist[[10]],scRNAlist[[11]],scRNAlist[[12]],scRNAlist[[13]],
                                           scRNAlist[[14]],scRNAlist[[15]],scRNAlist[[16]]))
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:25)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:25) %>% FindClusters()

plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T,raster=FALSE) 
plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident',raster=FALSE) 
plotc <- plot1+plot2
plotc
