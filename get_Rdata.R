## This file is used to convert the matrix into Rdata through seurat.

path = "~/Documents/Data/"
library(Seurat)
library(tidyverse)
E_MTAB_3929 <- read.delim(str_c(path, "scRecover/E-MTAB-3929/counts.txt"), row.names=1)
save(E_MTAB_3929, file = str_c(path, "scRecover/E-MTAB-3929/E_MTAB_3929.Rdata"))
E_MTAB_3929 <- CreateSeuratObject(counts = counts)

E_MTAB_3929 <- NormalizeData(object = E_MTAB_3929)
E_MTAB_3929 <- FindVariableFeatures(object = E_MTAB_3929)
E_MTAB_3929 <- ScaleData(object = E_MTAB_3929)
E_MTAB_3929 <- RunPCA(object = E_MTAB_3929)
E_MTAB_3929 <- FindNeighbors(object = E_MTAB_3929, dims = 1:30)
E_MTAB_3929 <- FindClusters(object = E_MTAB_3929)
#E_MTAB_3929 <- RunUMAP(object = E_MTAB_3929, dims = 1:30)
E_MTAB_3929 <- RunTSNE(object = E_MTAB_3929, dims = 1:30)
DimPlot(object = E_MTAB_3929, reduction = "tsne")

#save(E_MTAB_3929, file = str_c(path, "scRecover/E-MTAB-3929/E_MTAB_3929.Rdata"))
