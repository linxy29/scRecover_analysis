library(Seurat)
library(ggplot2)
library(clusteval)
library(mclust)
library(tidyverse)

seuratPipeline <- function(plotFolder, raw_data, dataName, label, label_2 = NULL, label_2_name = NULL, pcs.use = 20, resolution_min = 0.3, resolution_max = 1){
  #raw_data[raw_data < 0] <- 0
  S.raw_data <- CreateSeuratObject(counts =  raw_data)
  S.raw_data <- NormalizeData(object = S.raw_data)
  S.raw_data <- FindVariableFeatures(S.raw_data, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(S.raw_data), 10)
  plot1 <- VariableFeaturePlot(S.raw_data)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  ggsave(file = paste0(plotFolder, "VariableFeaturePlot of ", dataName, ".pdf"))
  all.genes <- rownames(S.raw_data)
  S.raw_data <- ScaleData(S.raw_data, features = all.genes)
  S.raw_data <- RunPCA(S.raw_data, features = VariableFeatures(object = S.raw_data))
  #DimHeatmap(S.raw_data, dims = 1:15, cells = 500, balanced = TRUE)
  #ggsave(file = paste0(plotFolder, "Top genes associated with principal components of ", dataName, ".pdf"))
  DimPlot(S.raw_data, reduction = "pca") + NoLegend()
  ggsave(file = paste0(plotFolder, "PCA of ", dataName, ".pdf"))
  ElbowPlot(S.raw_data)
  ggsave(file = paste0(plotFolder, "ElbowPlot of ", dataName, ".pdf"))
  
  
  # Plot t-SNE
  # while(length(unique(label[,1])) != length(levels(S.raw_data@ident))){
  #   cat("\r", "resolution = ", resolution)
  #   if(length(levels(S.raw_data@ident)) > length(unique(label[,1])))
  #     resolution <- resolution - 0.002
  #   if(length(levels(S.raw_data@ident)) < length(unique(label[,1])))
  #     resolution <- resolution + 0.05
  #   S.raw_data <- FindClusters(object = S.raw_data, reduction.type = "pca", dims.use = 1:pcs.use, resolution = resolution, print.output = FALSE, save.SNN = TRUE, force.recalc = TRUE)
  # }
  
  for(resolution in seq(resolution_min, resolution_max, by = 0.1)){
    S.raw_data <- FindNeighbors(S.raw_data, dims = 1:pcs.use)
    S.raw_data <- FindClusters(S.raw_data, resolution = resolution)
    S.raw_data <- RunTSNE(S.raw_data, dims = 1:pcs.use)
    
    seuratClust.ident <- Idents(S.raw_data)
    #saveRDS(seuratClust.ident, file = paste0(plotFolder, "seuratClust_", dataName, ".rds"))
    #plotConfusionMat(label, seuratClust.ident, plotFolder, dataName)
    
    DimPlot(S.raw_data, reduction = "tsne")
    ggsave(file = paste0(plotFolder, "t-SNE of ", dataName, " by seurat clusters_resolution_", resolution, ".pdf"), width = 6, height = 5)
    
    jaccard <- cluster_similarity(as.factor(label), seuratClust.ident, similarity="jaccard")
    indices <-  adjustedRandIndex(as.numeric(seuratClust.ident), as.numeric(label))
    methodName <- switch(EXPR = dataName, raw_data = "Raw data", scImpute_filter = "scRecover + scImpute", SAVER_filter = "scRecover + SAVER", MAGIC_filter = "scRecover + MAGIC", scimpute_count = "scImpute", SAVER_count = "SAVER", MAGIC_count = "MAGIC")
    indices_title <- paste(c("ARI", "Jaccard"), round(c(indices, jaccard), digits = 3), sep = " = ", collapse = ", ")
    indices_title <- paste0(methodName, "\n", indices_title)
    info <- paste0(indices_title, ", resolution = ", resolution, ", nc = ", length(levels(Idents(S.raw_data))))
    print(info)
    write(info, file = paste0(plotFolder, "info_", dataName, ".txt"), append = TRUE)
    Idents(S.raw_data) <- label
    
    DimPlot(S.raw_data, reduction = "tsne")  +
      ggtitle(paste0(methodName, "\nARI = ", as.character(round(indices, 3)), " JACCARD = ", as.character(round(jaccard, 3))))
    ggsave(file = paste0(plotFolder, "t-SNE of ", dataName, " by original clusters_resolution_", resolution, ".pdf"), width = 6, height = 5)
    
    if(!is.null(label_2)){
      jaccard <- cluster_similarity(as.factor(label_2[,1]), seuratClust.ident, similarity="jaccard")
      indices <- adjustedRand(as.integer(label_2[,1]), as.vector(seuratClust.ident), randMethod = c("HA", "Jaccard"))
      methodName <- switch(EXPR = dataName, raw_data = "Raw data", scImpute_filter = "scRecover + scImpute", SAVER_filter = "scRecover + SAVER", MAGIC_filter = "scRecover + MAGIC", scimpute_count = "scImpute", SAVER_count = "SAVER", MAGIC_count = "MAGIC")
      indices_title <- paste(c("ARI", "Jaccard"), round(indices, digits = 3), sep = " = ", collapse = ", ")
      indices_title <- paste0(methodName, "\n", indices_title)
      info <- paste0(indices_title, ", resolution = ", resolution, ", nc = ", length(levels(S.raw_data@ident)))
      write(info, file = paste0(plotFolder, "info_label_2_", dataName, ".txt"), append = TRUE)
      
      new.ident <- as.factor(label_2[,1])
      names(new.ident) <- row.names(label_2)
      S.raw_data@ident <- new.ident[colnames(S.raw_data@data)]
      
      png(paste0(plotFolder, "t-SNE of ", dataName, " by original clusters_", label_2_name, "_resolution_", resolution, ".png"), width = 325, height = 300)
      TSNEPlot(object = S.raw_data, pt.size = 2, do.label = TRUE, label.size = 6, plot.title = indices_title)
      dev.off()
    }
  }
  markers <- FindAllMarkers(S.raw_data, only.pos = TRUE)
  write_csv(markers, file = paste0(plotFolder, "markers of ", dataName, " by original clusters.csv"))
}

imputationAnalysis <- function(path, subplotFolder, label, label_2 = NULL, label_2_name = NULL, pcs.use = 20, resolution_min, resolution_max){
  plotFolder <- paste0(path, subplotFolder)
  dir.create(plotFolder, showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(plotFolder, "tempFile/"), showWarnings = FALSE, recursive = TRUE)
  ## read data
  counts_down <- read.csv(file = paste0(path, "/raw_data.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  #whether_impute_inz <- readRDS(file = paste0(path, "/tempFile/whether_impute_inz.rds"))
  scImpute_count <- read.csv(file = paste0(path, "/tempFile/scimpute_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  SAVER_count <- read.csv(file = paste0(path, "/tempFile/SAVER_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  #MAGIC_count <- read.csv(file = paste0(path, "/tempFile/MAGIC_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  scImpute_filter <- read.csv(file = paste0(path, "/scRecover+scImpute.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  SAVER_filter <- read.csv(file = paste0(path, "/scRecover+SAVER.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  #MAGIC_filter <- read.csv(file = paste0(path, "/scRecover+MAGIC.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  
  seuratPipeline(plotFolder = paste0(plotFolder, "tempFile/"), raw_data = counts_down, dataName = "counts_down", label = label, label_2 = label_2, label_2_name = label_2_name, pcs.use = pcs.use, resolution_min = resolution_min, resolution_max = resolution_max)
  seuratPipeline(plotFolder = plotFolder, raw_data = scImpute_filter, dataName = "scImpute_filter", label = label, label_2 = label_2, label_2_name = label_2_name, pcs.use = pcs.use, resolution_min = resolution_min, resolution_max = resolution_max)
  seuratPipeline(plotFolder = plotFolder, raw_data = SAVER_filter, dataName = "SAVER_filter", label = label, label_2 = label_2, label_2_name = label_2_name, pcs.use = pcs.use, resolution_min = resolution_min, resolution_max = resolution_max)
  #seuratPipeline(plotFolder = plotFolder, dataName = "MAGIC_filter", label = label, label_2 = label_2, label_2_name = label_2_name, pcs.use = pcs.use, resolution_max = resolution_max)
  seuratPipeline(plotFolder = paste0(plotFolder, "tempFile/"), raw_data = scImpute_count, dataName = "scimpute_count", label = label, label_2 = label_2, label_2_name = label_2_name, pcs.use = pcs.use, resolution_min = resolution_min, resolution_max = resolution_max)
  seuratPipeline(plotFolder = paste0(plotFolder, "tempFile/"), raw_data = SAVER_count, dataName = "SAVER_count", label = label, label_2 = label_2, label_2_name = label_2_name, pcs.use = pcs.use, resolution_min = resolution_min, resolution_max = resolution_max)
  #seuratPipeline(plotFolder = paste0(plotFolder, "tempFile/"), dataName = "MAGIC_count", label = label, label_2 = label_2, label_2_name = label_2_name, pcs.use = pcs.use, resolution_min = resolution_min, resolution_max = resolution_max)
}

# pbmc3k
setwd("~/Documents/Data/scRecover/pbmc3k")
pbmc = readRDS("./pbmc3k_seuratLable.rds")
label = Idents(pbmc)
# get the plot of original data
counts_original <- pbmc@assays$RNA@counts
raw_data <- as.matrix(counts_original)
dir.create("./seurat/", showWarnings = FALSE, recursive = TRUE)
seuratPipeline(plotFolder = "./seurat/", raw_data = raw_data, dataName = "raw_data", label = label, pcs.use = 20, resolution_min = 0.3, resolution_max = 0.3)

# get the plot of imputation data
for(Kcluster in c(2, 5, 10)){
  for(depth in c(10, 20, 50)){
    path <- paste0("./pbmc3k_K_", Kcluster, "_d_", depth, "/")
    print(paste0("|+++++++++++++++++++++Processing ", path))
    imputationAnalysis(path = path, subplotFolder = "seurat/", label = label, resolution_min = 0.1, resolution_max = 1)
    gc()
  }
}

######## Check markers
# Load packages
library(tidyverse)
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("~/Documents/Data/scRecover/pbmc3k")
scImpute_markers = read_csv("~/Documents/Data/scRecover/pbmc3k/pbmc3k_K_2_d_10/seurat/markers of scImpute_filter by original clusters.csv")
SAVER_markers = read_csv("~/Documents/Data/scRecover/pbmc3k/pbmc3k_K_2_d_10/seurat/markers of SAVER_filter by original clusters.csv")
original_markers = read_csv("~/Documents/Data/scRecover/pbmc3k/seurat/markers of raw_data by original clusters_resolution_0.3.csv")

perform_GO_analysis <- function(gene_sets, index_to_analyze, OrgDb, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE) {
  # Ensure index_to_analyze is within bounds
  if (index_to_analyze > length(gene_sets) || index_to_analyze < 1) {
    stop("Index to analyze is out of bounds.")
  }
  
  # Compute the unique genes in the specified gene set
  other_indices <- setdiff(1:length(gene_sets), index_to_analyze)
  all_other_genes <- unique(unlist(gene_sets[other_indices]))
  unique_genes <- setdiff(gene_sets[[index_to_analyze]], all_other_genes)
  
  if (length(unique_genes) > 0) {
    ego <- enrichGO(
      gene = unique_genes,
      OrgDb = OrgDb,
      keyType = keyType,
      ont = ont,
      pAdjustMethod = pAdjustMethod,
      qvalueCutoff = qvalueCutoff,
      readable = readable
    )
    
    if (length(ego) > 0) {
      print(summary(ego))
    } else {
      print(paste("No significant GO terms found for unique genes in dataset", index_to_analyze))
    }
  } else {
    print(paste("No unique genes in dataset", index_to_analyze))
  }
}

analyze_cluster_genes <- function(df1, df2, df3, cluster_name) {
  # Extract genes for the specified cluster from each dataframe
  genes_df1 <- filter(df1, cluster == cluster_name) %>% pull(gene)
  genes_df2 <- filter(df2, cluster == cluster_name) %>% pull(gene)
  genes_df3 <- filter(df3, cluster == cluster_name) %>% pull(gene)
  
  dir_path <- "./marker_analysis"
  
  # Check if the directory exists; if not, create it
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  venn.plot <- venn.diagram(
    x = list(DF1 = genes_df1, DF2 = genes_df2, DF3 = genes_df3),
    category.names = c("original", "SAVER + scRecover", "scImpute + scRecover"),
    filename = NULL
  )
  grid.draw(venn.plot)
  
  # Perform GO term analysis 
  gene_sets <- list(genes_df1, genes_df2, genes_df3)
  saveRDS(gene_sets, file = paste0("./marker_analysis/geneset_", cluster_name, ".rds"))
  print("Enriched pathway specific to original dataset:")
  ego1 = perform_GO_analysis(gene_sets, 1, org.Hs.eg.db)
  print("Enriched pathway specific to SAVER + scRecover dataset:")
  ego2 = perform_GO_analysis(gene_sets, 2, org.Hs.eg.db)
  print("Enriched pathway specific to scImpute + scRecover original dataset:")
  ego3 = perform_GO_analysis(gene_sets, 3, org.Hs.eg.db)
  ego_df <- rbind(as.data.frame(ego1) %>% mutate(source = "origin specific"), as.data.frame(ego2) %>% mutate(source = "SAVER + scRecover specific"), as.data.frame(ego3) %>% mutate(source = "scImpute + scRecover specific"))
  write.csv(ego_df, paste0("./marker_analysis/goterm_", cluster_name, ".csv"), row.names = FALSE)
  
}

for (celltype in unique(original_markers$cluster)) {
  analyze_cluster_genes(original_markers, SAVER_markers, scImpute_markers, celltype)
}



###################
setwd("~/Documents/Data/scRecover/pbmc3k")
pbmc = readRDS("./pbmc3k_seuratLable.rds")
label = Idents(pbmc)
expression_original <- pbmc@assays$RNA@counts
seurat_origin = CreateSeuratObject(counts =  expression_original, project = "origin")
seurat_origin$celltype = label
expression_SAVER = read.csv("~/Documents/Data/scRecover/pbmc3k/pbmc3k_K_2_d_10/scRecover+SAVER.csv", stringsAsFactors = FALSE, header = T, row.names = 1)
seurat_SAVER = CreateSeuratObject(counts =  expression_SAVER, project = "scRecover_SAVER")
seurat_SAVER$celltype = label
expression_scImpute = read.csv("~/Documents/Data/scRecover/pbmc3k/pbmc3k_K_2_d_10/scRecover+scImpute.csv", stringsAsFactors = FALSE, header = T, row.names = 1)
seurat_scImpute = CreateSeuratObject(counts =  expression_scImpute, project = "scRecover_scImpute")
seurat_scImpute$celltype = label

pbmc.combined <- merge(seurat_origin, y = c(seurat_SAVER, seurat_scImpute), add.cell.ids = c("origin", "scRecover_SAVER", "scRecover_scImpute"))

## check marker
pbmc.subset <- subset(pbmc.combined, subset = celltype == "Naive CD4 T")
geneset = readRDS("./marker_analysis/geneset_Naive CD4 T.rds")
diff1 = setdiff(geneset[[1]], c(geneset[[2]], geneset[[3]]))
diff2 = setdiff(geneset[[2]], c(geneset[[1]], geneset[[3]]))
diff3 = setdiff(geneset[[3]], c(geneset[[2]], geneset[[1]]))
Idents(pbmc.subset) = pbmc.subset$orig.ident
markers = FindAllMarkers(pbmc.subset, only.pos = TRUE)
#plot_gene = intersect(c(diff1, diff2, diff3), markers$gene)
#plot_gene = c(plot_gene, "IMPDH2", "ADSL")
#plot_gene = c(plot_gene, markers$gene[1:7])
plot_gene = c(diff1[2:3], "RELB", "TNFAIP8", "SMC2", "RPL29", "PIM1", "IL32", "PHF1", "GATA3", "IMPDH2", "ADSL")

DotPlot(pbmc.subset, features = plot_gene)
ggsave("figure6B.pdf", width = 12, height = 4)


