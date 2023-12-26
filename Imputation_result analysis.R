# HCA_CM_TPM
library(Seurat)
library(clues)
library(clusteval)

convertFile <- function(path){
  raw_data <- read.csv(file = paste0(path, "raw_data.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  saveRDS(raw_data, file = paste0(path, "raw_data.rds"))
  scImpute_filter <- read.csv(file = paste0(path, "scImpute_filter.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  saveRDS(scImpute_filter, file = paste0(path, "scImpute_filter.rds"))
  SAVER_filter <- read.csv(file = paste0(path, "SAVER_filter.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  saveRDS(SAVER_filter, file = paste0(path, "SAVER_filter.rds"))
  MAGIC_filter <- read.csv(file = paste0(path, "MAGIC_filter.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  saveRDS(MAGIC_filter, file = paste0(path, "MAGIC_filter.rds"))
  scimpute_count <- read.csv(file = paste0(path, "tempFile/scimpute_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  saveRDS(scimpute_count, file = paste0(path, "tempFile/scimpute_count.rds"))
  SAVER_count <- read.csv(file = paste0(path, "tempFile/SAVER_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  saveRDS(SAVER_count, file = paste0(path, "tempFile/SAVER_count.rds"))
  MAGIC_count <- read.csv(file = paste0(path, "tempFile/MAGIC_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  saveRDS(MAGIC_count, file = paste0(path, "tempFile/MAGIC_count.rds"))
}

# PLot confusion matrix function
plotConfusionMat <- function(label, seuratClust, plotFolder, dataName){
  library(ggplot2)
  confusionMat <- table(label[,1], seuratClust)
  confusionMat <- confusionMat[,order(max.col(t(confusionMat), 'first'))]
  confusionMat_log <- log(confusionMat[nrow(confusionMat):1,] + 1)
  methodName <- switch(EXPR = dataName, raw_data = "Raw data", scImpute_filter = "scRecover + scImpute", SAVER_filter = "scRecover + SAVER", MAGIC_filter = "scRecover + MAGIC", scimpute_count = "scImpute", SAVER_count = "SAVER", MAGIC_count = "MAGIC")
  
  ggplot(data = as.data.frame(confusionMat_log), aes(seuratClust, Var1, fill = Freq))+
    ggtitle(methodName)+
    geom_tile(color = "gray")+
    scale_fill_gradient2(low = "white", high = "green", mid = "red", midpoint = max(confusionMat_log)/2, limit = c(0,max(confusionMat_log)), space = "Lab", name="Color bar") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 0, vjust = 2, size = 10, hjust = 0.5),
          axis.text.y = element_text(vjust = 0.5, size = 10, hjust = 1))+
    coord_fixed()+
    geom_text(data=as.data.frame(confusionMat[nrow(confusionMat):1,]), aes(seuratClust, Var1, label = Freq), color = "black", size = 3) +
    theme(
      plot.title = element_text(size=10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.justification = "center",
      legend.position = "none",
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 4, barheight = 0.8, title.position = "top", title.hjust = 0.5))
  ggsave(paste0(plotFolder, "ConfusionMat_", dataName, ".pdf"), width = 3, height = 3)
}

seuratPipeline <- function(path, plotFolder, dataName, label, label_2 = NULL, label_2_name = NULL, pcs.use = 20, resolution_max = 1){
  raw_data <- readRDS(file = paste0(path, dataName, ".rds"))
  raw_data[raw_data < 0] <- 0
  
  S.raw_data <- CreateSeuratObject(raw.data = raw_data, project = dataName)
  S.raw_data <- NormalizeData(object = S.raw_data, normalization.method = "LogNormalize", scale.factor = 10000)
  png(paste0(plotFolder, "Variable genes of ", dataName, ".png"), width = 800, height = 600)
  S.raw_data <- FindVariableGenes(object = S.raw_data, mean.function = ExpMean, dispersion.function = LogVMR)
  dev.off()
  length(x = S.raw_data@var.genes)
  S.raw_data <- ScaleData(object = S.raw_data, vars.to.regress = c("nUMI"), do.par = TRUE, num.cores = 7)
  
  S.raw_data <- RunPCA(object = S.raw_data, pc.genes = S.raw_data@var.genes, do.print = FALSE, pcs.print = 1:5, genes.print = 5, pcs.compute = 25)
  png(paste0(plotFolder, "Top genes associated with principal components of ", dataName, ".png"), width = 800, height = 600)
  VizPCA(object = S.raw_data, pcs.use = 1:2, font.size = 0.8)
  dev.off()
  png(paste0(plotFolder, "PCA of ", dataName, ".png"), width = 600, height = 400)
  PCAPlot(object = S.raw_data, dim.1 = 1, dim.2 = 2, pt.size = 2)
  dev.off()
  S.raw_data <- ProjectPCA(object = S.raw_data, do.print = FALSE)
  png(paste0(plotFolder, "PC heatmap of ", dataName, ".png"), width = 1000, height = 1200)
  PCHeatmap(object = S.raw_data, pc.use = 1:25, cells.use = 100, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE, cexRow = 1.2)
  dev.off()
  S.raw_data <- JackStraw(object = S.raw_data, num.pc = 25, num.replicate = 100, prop.freq = 0.01, do.par = TRUE, num.cores = 7)
  png(paste0(plotFolder, "Jack Straw Plot of ", dataName, ".png"), width = 700, height = 700)
  JackStrawPlot(object = S.raw_data, PCs = 1:25, nCol = 5)
  dev.off()
  png(paste0(plotFolder, "Standard deviations of PCs of ", dataName, ".png"), width = 400, height = 400)
  PCElbowPlot(object = S.raw_data, num.pc = 25)
  dev.off()
  
  # Plot t-SNE
  # while(length(unique(label[,1])) != length(levels(S.raw_data@ident))){
  #   cat("\r", "resolution = ", resolution)
  #   if(length(levels(S.raw_data@ident)) > length(unique(label[,1])))
  #     resolution <- resolution - 0.002
  #   if(length(levels(S.raw_data@ident)) < length(unique(label[,1])))
  #     resolution <- resolution + 0.05
  #   S.raw_data <- FindClusters(object = S.raw_data, reduction.type = "pca", dims.use = 1:pcs.use, resolution = resolution, print.output = FALSE, save.SNN = TRUE, force.recalc = TRUE)
  # }
  
  for(resolution in seq(resolution_max, resolution_max, by = 0.1)){
    
    S.raw_data <- FindClusters(object = S.raw_data, reduction.type = "pca", dims.use = 1:pcs.use, resolution = resolution, print.output = FALSE, save.SNN = TRUE, force.recalc = TRUE)
    S.raw_data <- RunTSNE(object = S.raw_data, dims.use = 1:pcs.use, do.fast = TRUE, perplexity = 30, seed.use = 1)
    seuratClust.ident <- S.raw_data@ident
    saveRDS(seuratClust.ident, file = paste0(plotFolder, "seuratClust_", dataName, ".rds"))
    plotConfusionMat(label, seuratClust.ident, plotFolder, dataName)
    
    png(paste0(plotFolder, "t-SNE of ", dataName, " by seurat clusters_resolution_", resolution, ".png"), width = 325, height = 300)
    TSNEPlot(object = S.raw_data, pt.size = 2, do.label = TRUE, label.size = 6)
    dev.off()
    
    jaccard <- cluster_similarity(as.factor(label[,1]), seuratClust.ident, similarity="jaccard")
    indices <- adjustedRand(as.integer(label[,1]), as.vector(seuratClust.ident), randMethod = c("HA", "Jaccard"))
    methodName <- switch(EXPR = dataName, raw_data = "Raw data", scImpute_filter = "scRecover + scImpute", SAVER_filter = "scRecover + SAVER", MAGIC_filter = "scRecover + MAGIC", scimpute_count = "scImpute", SAVER_count = "SAVER", MAGIC_count = "MAGIC")
    indices_title <- paste(c("ARI", "Jaccard"), round(indices, digits = 3), sep = " = ", collapse = ", ")
    indices_title <- paste0(methodName, "\n", indices_title)
    info <- paste0(indices_title, ", resolution = ", resolution, ", nc = ", length(levels(S.raw_data@ident)))
    write(info, file = paste0(plotFolder, "info_", dataName, ".txt"), append = TRUE)
    
    new.ident <- as.factor(label[,1])
    names(new.ident) <- row.names(label)
    S.raw_data@ident <- new.ident[colnames(S.raw_data@data)]
    
    png(paste0(plotFolder, "t-SNE of ", dataName, " by original clusters_resolution_", resolution, ".png"), width = 325, height = 300)
    TSNEPlot(object = S.raw_data, pt.size = 2, do.label = TRUE, label.size = 6, plot.title = indices_title)
    dev.off()
    
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
}

imputationAnalysis <- function(path, plotFolder, label, label_2 = NULL, label_2_name = NULL, pcs.use = 20, resolution_max = 1){
  plotFolder <- paste0(path, plotFolder)
  dir.create(plotFolder, showWarnings = FALSE, recursive = TRUE)
  seuratPipeline(path = path, plotFolder = plotFolder, dataName = "raw_data", label = label, label_2 = label_2, label_2_name = label_2_name, pcs.use = pcs.use, resolution_max = resolution_max)
  seuratPipeline(path = path, plotFolder = plotFolder, dataName = "scImpute_filter", label = label, label_2 = label_2, label_2_name = label_2_name, pcs.use = pcs.use, resolution_max = resolution_max)
  seuratPipeline(path = path, plotFolder = plotFolder, dataName = "SAVER_filter", label = label, label_2 = label_2, label_2_name = label_2_name, pcs.use = pcs.use, resolution_max = resolution_max)
  seuratPipeline(path = path, plotFolder = plotFolder, dataName = "MAGIC_filter", label = label, label_2 = label_2, label_2_name = label_2_name, pcs.use = pcs.use, resolution_max = resolution_max)
  seuratPipeline(path = paste0(path, "tempFile/"), plotFolder = plotFolder, dataName = "scimpute_count", label = label, label_2 = label_2, label_2_name = label_2_name, pcs.use = pcs.use, resolution_max = resolution_max)
  seuratPipeline(path = paste0(path, "tempFile/"), plotFolder = plotFolder, dataName = "SAVER_count", label = label, label_2 = label_2, label_2_name = label_2_name, pcs.use = pcs.use, resolution_max = resolution_max)
  seuratPipeline(path = paste0(path, "tempFile/"), plotFolder = plotFolder, dataName = "MAGIC_count", label = label, label_2 = label_2, label_2_name = label_2_name, pcs.use = pcs.use, resolution_max = resolution_max)
}



library(ggplot2)
library(RColorBrewer)
imputationBasicPlot <- function(path){
  plotFolder <- paste0(path, "Basic Plot/")
  dir.create(plotFolder, showWarnings = FALSE)
  raw_data <- readRDS(file = paste0(path, "raw_data.rds"))
  whether_impute_inz <- readRDS(file = paste0(path, "tempFile/whether_impute_inz.rds"))
  scImpute_count <- readRDS(file = paste0(path, "tempFile/scimpute_count.rds"))
  SAVER_count <- readRDS(file = paste0(path, "tempFile/SAVER_count.rds"))
  MAGIC_count <- readRDS(file = paste0(path, "tempFile/MAGIC_count.rds"))
  scImpute_filter <- readRDS(file = paste0(path, "scImpute_filter.rds"))
  SAVER_filter <- readRDS(file = paste0(path, "SAVER_filter.rds"))
  MAGIC_filter <- readRDS(file = paste0(path, "MAGIC_filter.rds"))
  
  # Gene number barplot
  countGenes <- function(raw_data, counts){
    return(c(mean(colSums(raw_data != 0)), mean(colSums(raw_data == 0 & counts != 0)), mean(colSums(raw_data == 0 & counts == 0))))
  }
  geneNumberMatrix <- cbind(countGenes(raw_data, raw_data), countGenes(raw_data, SAVER_count), countGenes(raw_data, SAVER_filter), countGenes(raw_data, scImpute_count), countGenes(raw_data, scImpute_filter), countGenes(raw_data, MAGIC_count), countGenes(raw_data, MAGIC_filter))
  row.names(geneNumberMatrix) <- c("Expressed", "Imputed", "Zeros")
  colnames(geneNumberMatrix) <- c("Raw data", "SAVER", "scRecover + SAVER", "scImpute", "scRecover + scImpute", "MAGIC", "scRecover + MAGIC")
  png(paste0(plotFolder, "Gene numbers.png"), width = 500, height = 500)
  col = brewer.pal(3, "Pastel2")
  b <- barplot(geneNumberMatrix, col=col , border="white", ylim = c(0, 1000*ceiling(colSums(geneNumberMatrix)[1]/1000)), xaxt="n")
  text(cex=1.2, x= b + 0.1, y = par("usr")[3], labels = colnames(geneNumberMatrix), srt=30, adj=1, xpd=TRUE)
  labelText <- as.character(round(as.vector(t(geneNumberMatrix))))
  labelText[labelText == "0"] <- NA
  text(cex=1.2, c(b, b, b), c(geneNumberMatrix[1,]/2, geneNumberMatrix[1,] + geneNumberMatrix[2,]/2, geneNumberMatrix[1,] + geneNumberMatrix[2,] + geneNumberMatrix[3,]/2), labelText)
  par(xpd=TRUE)
  legend(0.5, 1000*ceiling(colSums(geneNumberMatrix)[1]/1000)*1.08, legend = row.names(geneNumberMatrix), fill = col, horiz = TRUE, bty = "n", border = "white", cex = 1.4)
  dev.off()
  
  # Gene number distribution
  GN_raw_data <- colSums(raw_data != 0)
  GN_whether_impute <- colSums(whether_impute_inz != 0)
  GN_scImpute_count <- colSums(scImpute_count != 0)
  GN_scImpute_filter <- colSums(scImpute_filter != 0)
  GN_SAVER_count <- colSums(SAVER_count != 0)
  GN_SAVER_filter <- colSums(SAVER_filter != 0)
  GN_MAGIC_count <- colSums(MAGIC_count != 0)
  GN_MAGIC_filter <- colSums(MAGIC_filter != 0)
  GN_all <- data.frame(GeneNumber=c(GN_raw_data, GN_whether_impute, GN_scImpute_count, GN_scImpute_filter, GN_SAVER_count, GN_SAVER_filter, GN_MAGIC_count, GN_MAGIC_filter), method=rep(c("Raw", "scRecover", "scImpute", "scImpute+", "SAVER", "SAVER+", "MAGIC", "MAGIC+"), each = ncol(raw_data)))
  ggplot(data=GN_all, aes(x=GeneNumber, colour=method, group=method)) + geom_density(alpha=0.5, size = 1.2)
  ggsave(paste0(plotFolder, "Gene number distribution.png"), width = 5, height = 4)
  
  GN_raw_data <- rowMeans(raw_data != 0)
  GN_whether_impute <- rowMeans(whether_impute_inz != 0)
  GN_scImpute_count <- rowMeans(scImpute_count != 0)
  GN_scImpute_filter <- rowMeans(scImpute_filter != 0)
  GN_SAVER_count <- rowMeans(SAVER_count != 0)
  GN_SAVER_filter <- rowMeans(SAVER_filter != 0)
  GN_MAGIC_count <- rowMeans(MAGIC_count != 0)
  GN_MAGIC_filter <- rowMeans(MAGIC_filter != 0)
  GN_all <- data.frame(ExpressedRatio=c(GN_raw_data, GN_whether_impute, GN_scImpute_count, GN_scImpute_filter, GN_SAVER_count, GN_SAVER_filter, GN_MAGIC_count, GN_MAGIC_filter), method=rep(c("Raw", "scRecover", "scImpute", "scImpute+", "SAVER", "SAVER+", "MAGIC", "MAGIC+"), each = nrow(raw_data)))
  ggplot(data=GN_all, aes(x=ExpressedRatio, colour=method, group=method)) + geom_density(alpha=0.5, size = 1.2)
  ggsave(paste0(plotFolder, "Gene expressed ratio distribution.png"), width = 5, height = 4)
  
  GN_raw_data <- colSums(raw_data);GN_raw_data <- GN_raw_data[GN_raw_data < quantile(GN_raw_data, probs = 0.95)]
  GN_scImpute_count <- colSums(scImpute_count);GN_scImpute_count <- GN_scImpute_count[GN_scImpute_count < quantile(GN_scImpute_count, probs = 0.95)]
  GN_scImpute_filter <- colSums(scImpute_filter);GN_scImpute_filter <- GN_scImpute_filter[GN_scImpute_filter < quantile(GN_scImpute_filter, probs = 0.95)]
  GN_SAVER_count <- colSums(SAVER_count);GN_SAVER_count <- GN_SAVER_count[GN_SAVER_count < quantile(GN_SAVER_count, probs = 0.95)]
  GN_SAVER_filter <- colSums(SAVER_filter);GN_SAVER_filter <- GN_SAVER_filter[GN_SAVER_filter < quantile(GN_SAVER_filter, probs = 0.95)]
  GN_MAGIC_count <- colSums(MAGIC_count);GN_MAGIC_count <- GN_MAGIC_count[GN_MAGIC_count < quantile(GN_MAGIC_count, probs = 0.95)]
  GN_MAGIC_filter <- colSums(MAGIC_filter);GN_MAGIC_filter <- GN_MAGIC_filter[GN_MAGIC_filter < quantile(GN_MAGIC_filter, probs = 0.95)]
  GN_all <- data.frame(Librarysize=c(GN_raw_data, GN_scImpute_count, GN_scImpute_filter, GN_SAVER_count, GN_SAVER_filter, GN_MAGIC_count, GN_MAGIC_filter), method=rep(c("Raw", "scImpute", "scImpute+", "SAVER", "SAVER+", "MAGIC", "MAGIC+"), times = c(length(GN_raw_data), length(GN_scImpute_count), length(GN_scImpute_filter), length(GN_SAVER_count), length(GN_SAVER_filter), length(GN_MAGIC_count), length(GN_MAGIC_filter))))
  ggplot(data=GN_all, aes(x=Librarysize, colour=method, group=method)) + geom_density(alpha=0.5, size = 1.2)
  ggsave(paste0(plotFolder, "Library size distribution.png"), width = 5, height = 4)
  
  GN_raw_data <- rowMeans(raw_data);GN_raw_data <- GN_raw_data[GN_raw_data < quantile(GN_raw_data, probs = 0.95)]
  GN_scImpute_count <- rowMeans(scImpute_count);GN_scImpute_count <- GN_scImpute_count[GN_scImpute_count < quantile(GN_scImpute_count, probs = 0.95)]
  GN_scImpute_filter <- rowMeans(scImpute_filter);GN_scImpute_filter <- GN_scImpute_filter[GN_scImpute_filter < quantile(GN_scImpute_filter, probs = 0.95)]
  GN_SAVER_count <- rowMeans(SAVER_count);GN_SAVER_count <- GN_SAVER_count[GN_SAVER_count < quantile(GN_SAVER_count, probs = 0.95)]
  GN_SAVER_filter <- rowMeans(SAVER_filter);GN_SAVER_filter <- GN_SAVER_filter[GN_SAVER_filter < quantile(GN_SAVER_filter, probs = 0.95)]
  GN_MAGIC_count <- rowMeans(MAGIC_count);GN_MAGIC_count <- GN_MAGIC_count[GN_MAGIC_count < quantile(GN_MAGIC_count, probs = 0.95)]
  GN_MAGIC_filter <- rowMeans(MAGIC_filter);GN_MAGIC_filter <- GN_MAGIC_filter[GN_MAGIC_filter < quantile(GN_MAGIC_filter, probs = 0.95)]
  GN_all <- data.frame(MeanExpression=c(GN_raw_data, GN_scImpute_count, GN_scImpute_filter, GN_SAVER_count, GN_SAVER_filter, GN_MAGIC_count, GN_MAGIC_filter), method=rep(c("Raw", "scImpute", "scImpute+", "SAVER", "SAVER+", "MAGIC", "MAGIC+"), times = c(length(GN_raw_data), length(GN_scImpute_count), length(GN_scImpute_filter), length(GN_SAVER_count), length(GN_SAVER_filter), length(GN_MAGIC_count), length(GN_MAGIC_filter))))
  ggplot(data=GN_all, aes(x=MeanExpression, colour=method, group=method)) + geom_density(alpha=0.5, size = 1.2)
  ggsave(paste0(plotFolder, "Gene mean expression distribution.png"), width = 5, height = 4)
  
}








# 10x_neuron_1k_v3
label <- read.csv(file = "F:/data/10x_neuron_1k_v3/analysis/clustering/graphclust/clusters.csv", stringsAsFactors = FALSE, header = T, row.names = 1)
row.names(label) <- gsub("-", ".", row.names(label))

path <- "F:/Plots/ImputeSingle/10x_neuron_1k_v3/10x_neuron_1k_v3_ImputeSingle_full_K_5_d_20/"
imputationAnalysis(path = path, label = label)
imputationBasicPlot(path)

raw_data <- read.csv(file = paste0(path, "raw_data.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
whether_impute_inz <- readRDS(file = paste0(path, "tempFile/whether_impute_inz.rds"))
scImpute_count <- read.csv(file = paste0(path, "tempFile/scimpute_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
SAVER_count <- read.csv(file = paste0(path, "tempFile/SAVER_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
MAGIC_count <- read.csv(file = paste0(path, "tempFile/MAGIC_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
scImpute_filter <- read.csv(file = paste0(path, "scImpute_filter.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
SAVER_filter <- read.csv(file = paste0(path, "SAVER_filter.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
MAGIC_filter <- read.csv(file = paste0(path, "MAGIC_filter.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)

mean(colSums(whether_impute_inz != 0 & raw_data == 0))
mean(colSums(scImpute_count != 0 & raw_data == 0))
mean(colSums(SAVER_count != 0 & raw_data == 0))
mean(colSums(MAGIC_count != 0 & raw_data == 0))
mean(colSums(scImpute_filter != 0 & raw_data == 0))
mean(colSums(SAVER_filter != 0 & raw_data == 0))
mean(colSums(MAGIC_filter != 0 & raw_data == 0))


# 10x_pbmc_1k_v3
label <- read.csv(file = "F:/data/10x_pbmc_1k_v3/analysis/clustering/graphclust/clusters.csv", stringsAsFactors = FALSE, header = T, row.names = 1)
row.names(label) <- gsub("-", ".", row.names(label))

path <- "F:/Plots/ImputeSingle/10x_pbmc_1k_v3/10x_pbmc_1k_v3_ImputeSingle_full_K_5_d_20/"
imputationAnalysis(path = path, label = label, resolution = 0.5)
imputationBasicPlot(path)


# 10x_heart_1k_v3
label <- read.csv(file = "F:/data/10x_heart_1k_v3/analysis/clustering/graphclust/clusters.csv", stringsAsFactors = FALSE, header = T, row.names = 1)
row.names(label) <- gsub("-", ".", row.names(label))

path <- "F:/Plots/ImputeSingle/10x_heart_1k_v3/10x_heart_1k_v3_ImputeSingle_full_K_5_d_20/"
# convertFile(path)
# resolution = 0.3
imputationAnalysis(path = path, plotFolder = "tSNE/", label = label, resolution_max = 0.3)
imputationBasicPlot(path)


# E_MTAB_3929
load("F:/data/E-MTAB-3929/E_MTAB_3929.Rdata")
E_MTAB_3929 <- E_MTAB_3929[rowSums(E_MTAB_3929 != 0) > 0,]

SDRF_original <- read.table("F:/data/E-MTAB-3929/E-MTAB-3929.sdrf.txt", sep = "\t", header = T, stringsAsFactors = F)
SDRF <- SDRF_original[,c(1, 9, 10, 11, 12, 13, 14, 15)]
row.names(SDRF) <- SDRF[,1]
colnames(SDRF) <- c("sample.name", "individual", "developmental.stage", "treatment", "phenotype", "inferred.lineage", "inferred.TE.subpopulation", "inferred.pseudo.time")

label_day <- SDRF[match(colnames(E_MTAB_3929), row.names(SDRF)), "developmental.stage"]
label_lineage <- SDRF[match(colnames(E_MTAB_3929), row.names(SDRF)), "inferred.lineage"]
label_day[label_day == "embryonic day 3"] <- "E3"
label_day[label_day == "embryonic day 4"] <- "E4"
label_day[label_day == "embryonic day 5"] <- "E5"
label_day[label_day == "embryonic day 6"] <- "E6"
label_day[label_day == "embryonic day 7"] <- "E7"
label_lineage[label_lineage == "trophectoderm"] <- "TE"
label_lineage[label_lineage == "primitive endoderm"] <- "PE"
label_lineage[label_lineage == "epiblast"] <- "EPI"
label_lineage[label_lineage == "not applicable"] <- "NA"
label_day <- as.data.frame(label_day)
label_lineage <- as.data.frame(label_lineage)
row.names(label_day) <- colnames(E_MTAB_3929)
row.names(label_lineage) <- colnames(E_MTAB_3929)

path <- "F:/Plots/ImputeSingle/E_MTAB_3929/E_MTAB_3929_ImputeSingle_full_K_5_d_20/"
# convertFile(path)
imputationAnalysis(path = path, label = label_day, label_2 = label_lineage, label_2_name = "lineage", resolution = 0.1)
imputationBasicPlot(path)


# Deng
label <- readRDS(file = "F:/data/GSE45719/label_GSE45719.rds")
row.names(label) <- gsub("-", ".", row.names(label))
path <- "F:/Plots/ImputeSingle/Deng/Deng_ImputeSingle_full_K_5_d_20/"
# convertFile(path)
# resolution = 3.4
imputationAnalysis(path = path, plotFolder = "tSNE/", label = label, resolution_max = 5)
imputationBasicPlot(path)


# Chu1
label <- readRDS(file = "F:/data/GSE75748/GSE75748_sc_cell_type_label.rds")
path <- "F:/Plots/ImputeSingle/Chu1/Chu1_ImputeSingle_full_K_5_d_20/"
# convertFile(path)
# resolution = 0.5
imputationAnalysis(path = path, plotFolder = "tSNE/", label = label, resolution_max = 0.5)
imputationBasicPlot(path)


# Chu2
label <- readRDS(file = "F:/data/GSE75748/GSE75748_sc_time_course_label.rds")
path <- "F:/Plots/ImputeSingle/Chu2/Chu2_ImputeSingle_full_K_5_d_20/"
# convertFile(path)
imputationAnalysis(path = path, label = label, resolution = 0.5)
imputationBasicPlot(path)


# Li
label <- readRDS(file = "F:/data/GSE81861/GSE81861_Cell_Line_label.rds")
path <- "F:/Plots/ImputeSingle/Li/Li_ImputeSingle_full_K_5_d_20/"
# convertFile(path)
imputationAnalysis(path = path, label = label, resolution = 2.5)
imputationBasicPlot(path)









