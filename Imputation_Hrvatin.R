# ImputeSingle Hrvatin
library(ImputeSingle)
dat <- readRDS("F:/data/SAVER-data/hrvatin.rds")
set.seed(999)
dat <- dat[, sample(1:ncol(dat), 1000)]
ImputeSingle(counts = dat, Kcluster = 35, outputDir = "F:/Plots/ImputeSingle/Hrvatin_ImputeSingle_depth5/", depth = 5, SAVER = TRUE, MAGIC = TRUE, parallel = TRUE)

# Further analysis of Hrvatin
dat <- readRDS("F:/data/SAVER-data/hrvatin.rds")
set.seed(999)
dat <- dat[, sample(1:ncol(dat), 1000)]

load("F:/Plots/ImputeSingle/Hrvatin_ImputeSingle_depth5/tempFile/IntermediateVariables.Rdata")
x.magic <- read.csv(file = "F:/Plots/ImputeSingle/Hrvatin_ImputeSingle_depth5/tempFile/MAGIC_count.csv", stringsAsFactors = FALSE, header = T, row.names = 1)
x.saver <- read.csv(file = "F:/Plots/ImputeSingle/Hrvatin_ImputeSingle_depth5/tempFile/SAVER_count.csv", stringsAsFactors = FALSE, header = T, row.names = 1)
x.scImpute <- read.csv(file = "F:/Plots/ImputeSingle/Hrvatin_ImputeSingle_depth5/tempFile/scimpute_count.csv", stringsAsFactors = FALSE, header = T, row.names = 1)

x.magic_inz <- read.csv(file = "F:/Plots/ImputeSingle/Hrvatin_ImputeSingle_depth5/counts_MAGIC_inz.csv", stringsAsFactors = FALSE, header = T, row.names = 1)
x.saver_inz <- read.csv(file = "F:/Plots/ImputeSingle/Hrvatin_ImputeSingle_depth5/counts_SAVER_inz.csv", stringsAsFactors = FALSE, header = T, row.names = 1)
x.scImpute_inz <- read.csv(file = "F:/Plots/ImputeSingle/Hrvatin_ImputeSingle_depth5/counts_scImpute_inz.csv", stringsAsFactors = FALSE, header = T, row.names = 1)

# mean(colSums(dat == 0))
# mean(colSums(x.magic != 0 & dat == 0))
# mean(colSums(x.saver != 0 & dat == 0))
# mean(colSums(x.scImpute != 0 & dat == 0))
# mean(colSums(whether_impute_iz))
# 
# mean(colSums(x.magic != 0 & dat == 0 & whether_impute_iz == TRUE))
# mean(colSums(x.saver != 0 & dat == 0 & whether_impute_iz == TRUE))
# mean(colSums(x.scImpute != 0 & dat == 0 & whether_impute_iz == TRUE))

x <- dat

ct.dat <- read.csv("F:/data/SAVER-data/hrvatin_celltypes.csv", header = TRUE,
                   row.names = 1, stringsAsFactors = FALSE)

ct.dat2 <- ct.dat[!is.na(ct.dat$celltype), ]
ct.dat2$celltype[!is.na(ct.dat2$subtype)] <-
  ct.dat2$subtype[!is.na(ct.dat2$subtype)]
table(ct.dat2$celltype)
ct.dat3 <- ct.dat2[ct.dat2$celltype != "ExcL23", ]

ident <- ct.dat3[colnames(x), 4]

x.sub <- x[, !is.na(ident)]
x.saver.sub <- x.saver[, !is.na(ident)]
x.magic.sub <- x.magic[, !is.na(ident)]
x.scImpute.sub <- x.scImpute[, !is.na(ident)]
x.saver_inz.sub <- x.saver_inz[, !is.na(ident)]
x.magic_inz.sub <- x.magic_inz[, !is.na(ident)]
x.scImpute_inz.sub <- x.scImpute_inz[, !is.na(ident)]

ident2 <- ct.dat3[colnames(x.sub), 4]

non.exc <- grep("ExcL|Hip|RSP|Sub", ident2, invert = TRUE)

ident2[non.exc] <- "Other"

ident2 <- as.factor(ident2)

library(Seurat)
library(ggplot2)
library(cowplot)


obs <- CreateSeuratObject(x.sub, project = "obs")
saver <- CreateSeuratObject(x.saver.sub, project = "saver")
magic <- CreateSeuratObject(x.magic.sub, project = "magic")
scImpute <- CreateSeuratObject(x.scImpute.sub, project = "scImpute")
saver_inz <- CreateSeuratObject(x.saver_inz.sub, project = "saver_inz")
magic_inz <- CreateSeuratObject(x.magic_inz.sub, project = "magic_inz")
scImpute_inz <- CreateSeuratObject(x.scImpute_inz.sub, project = "scImpute_inz")

obs <- NormalizeData(obs)
saver <- NormalizeData(saver)
magic <- NormalizeData(magic)
scImpute <- NormalizeData(scImpute)
saver_inz <- NormalizeData(saver_inz)
magic_inz <- NormalizeData(magic_inz)
scImpute_inz <- NormalizeData(scImpute_inz)

obs <- FindVariableGenes(obs, do.plot = FALSE)
saver <- FindVariableGenes(saver, do.plot = FALSE)
magic <- FindVariableGenes(magic, do.plot = FALSE)
scImpute <- FindVariableGenes(scImpute, do.plot = FALSE)
saver_inz <- FindVariableGenes(saver_inz, do.plot = FALSE)
magic_inz <- FindVariableGenes(magic_inz, do.plot = FALSE)
scImpute_inz <- FindVariableGenes(scImpute_inz, do.plot = FALSE)

obs <- ScaleData(obs)
saver <- ScaleData(saver)
magic <- ScaleData(magic)
scImpute <- ScaleData(scImpute)
saver_inz <- ScaleData(saver_inz)
magic_inz <- ScaleData(magic_inz)
scImpute_inz <- ScaleData(scImpute_inz)

obs <- RunPCA(obs, pc.genes = obs@var.genes, do.print = FALSE, pcs.compute = 50)
saver <- RunPCA(saver, pc.genes = saver@var.genes, do.print = FALSE, pcs.compute = 50)
magic <- RunPCA(magic, pc.genes = magic@var.genes, do.print = FALSE, pcs.compute = 50)
scImpute <- RunPCA(scImpute, pc.genes = scImpute@var.genes, do.print = FALSE, pcs.compute = 50)
saver_inz <- RunPCA(saver_inz, pc.genes = saver_inz@var.genes, do.print = FALSE, pcs.compute = 50)
magic_inz <- RunPCA(magic_inz, pc.genes = magic_inz@var.genes, do.print = FALSE, pcs.compute = 50)
scImpute_inz <- RunPCA(scImpute_inz, pc.genes = scImpute_inz@var.genes, do.print = FALSE, pcs.compute = 50)

PCElbowPlot(obs, num.pc = 50)
PCElbowPlot(saver, num.pc = 50)
PCElbowPlot(magic, num.pc = 50)
PCElbowPlot(scImpute, num.pc = 50)
PCElbowPlot(saver_inz, num.pc = 50)
PCElbowPlot(magic_inz, num.pc = 50)
PCElbowPlot(scImpute_inz, num.pc = 50)

obs <- JackStraw(object = obs, num.pc = 40, num.replicate = 100)
saver <- JackStraw(object = saver, num.pc = 40, num.replicate = 100)
magic <- JackStraw(object = magic, num.pc = 40, num.replicate = 100)
scImpute <- JackStraw(object = scImpute, num.pc = 40, num.replicate = 100)
saver_inz <- JackStraw(object = saver_inz, num.pc = 40, num.replicate = 100)
magic_inz <- JackStraw(object = magic_inz, num.pc = 40, num.replicate = 100)
scImpute_inz <- JackStraw(object = scImpute_inz, num.pc = 40, num.replicate = 100)

p1 <- JackStrawPlot(obs, PCs = 1:40)
p2 <- JackStrawPlot(saver, PCs = 1:40)
p3 <- JackStrawPlot(magic, PCs = 1:40)
p4 <- JackStrawPlot(scImpute, PCs = 1:40)
p5 <- JackStrawPlot(saver_inz, PCs = 1:40)
p6 <- JackStrawPlot(magic_inz, PCs = 1:40)
p7 <- JackStrawPlot(scImpute_inz, PCs = 1:40)

levels(p1$data$PC.Score)
levels(p2$data$PC.Score)
levels(p3$data$PC.Score)
levels(p4$data$PC.Score)
levels(p5$data$PC.Score)
levels(p6$data$PC.Score)
levels(p7$data$PC.Score)

dim.obs <- 30
dim.saver <- 20
dim.magic <- 20
dim.scImpute <- 20
dim.saver_inz <- 20
dim.magic_inz <- 20
dim.scImpute_inz <- 20

obs <- RunTSNE(obs, dims.use = 1:dim.obs, check_duplicates = FALSE, do.fast = TRUE)
saver <- RunTSNE(saver, dims.use = 1:dim.saver, check_duplicates = FALSE, do.fast = TRUE)
magic <- RunTSNE(magic, dims.use = 1:dim.magic, check_duplicates = FALSE, do.fast = TRUE)
scImpute <- RunTSNE(scImpute, dims.use = 1:dim.scImpute, check_duplicates = FALSE, do.fast = TRUE)
saver_inz <- RunTSNE(saver_inz, dims.use = 1:dim.saver_inz, check_duplicates = FALSE, do.fast = TRUE)
magic_inz <- RunTSNE(magic_inz, dims.use = 1:dim.magic_inz, check_duplicates = FALSE, do.fast = TRUE)
scImpute_inz <- RunTSNE(scImpute_inz, dims.use = 1:dim.scImpute_inz, check_duplicates = FALSE, do.fast = TRUE)

ident3 <- factor(ident2, levels = levels(ident2)[c(3, 4, 5, 6, 7, 8, 9, 1, 2, 10, 12, 13, 11)])

obs.tsne <- data.frame(obs@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
saver.tsne <- data.frame(saver@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
magic.tsne <- data.frame(magic@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
scImpute.tsne <- data.frame(scImpute@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
saver_inz.tsne <- data.frame(saver_inz@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
magic_inz.tsne <- data.frame(magic_inz@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
scImpute_inz.tsne <- data.frame(scImpute_inz@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

my.cols <- c(gg_color_hue(12), "gray50")

g1 <- ggplot(obs.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2,
                                        colour = ident3), size = 0.5) +
  scale_colour_manual(name = "Cell type", values = my.cols) +
  guides(colour = FALSE) + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "plain")) +
  labs(title = "Observed")

g2 <- ggplot(saver.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2,
                                          colour = ident3), size = 0.5) +
  scale_colour_manual(name = "Cell type", values = my.cols) +
  guides(colour = FALSE) + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "plain")) +
  labs(title = "SAVER")

g3 <- ggplot(magic.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2,
                                          colour = ident3), size = 0.5) +
  scale_colour_manual(name = "Cell type", values = my.cols) +
  guides(colour = FALSE) + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "plain")) +
  labs(title = "MAGIC")

g4 <- ggplot(scImpute.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2,
                                             colour = ident3), size = 0.5) +
  scale_colour_manual(name = "Cell type", values = my.cols) +
  guides(colour = FALSE) + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "plain")) +
  labs(title = "scImpute")

g5 <- ggplot(saver_inz.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2,
                                              colour = ident3), size = 0.5) +
  scale_colour_manual(name = "Cell type", values = my.cols) +
  guides(colour = FALSE) + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "plain")) +
  labs(title = "SAVER_inz")

g6 <- ggplot(magic_inz.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2,
                                              colour = ident3), size = 0.5) +
  scale_colour_manual(name = "Cell type", values = my.cols) +
  guides(colour = FALSE) + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "plain")) +
  labs(title = "MAGIC_inz")

g7 <- ggplot(scImpute_inz.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2,
                                                 colour = ident3), size = 0.5) +
  scale_colour_manual(name = "Cell type", values = my.cols) +
  guides(colour = FALSE) + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "plain")) +
  labs(title = "scImpute_inz")


g0 <- ggplot(saver.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2,
                                          colour = ident3), size = 0.5) +
  scale_colour_manual(name = "Cell type", values = my.cols) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 14)) +
  labs(title = "SAVER")

gglegend <- function(x){ 
  tmp <- ggplot_gtable(ggplot_build(x)) 
  leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
  tmp$grobs[[leg]]
}

png("F:/Plots/ImputeSingle/fig2e_hrvatin.png", 9, 9, units = "in", res = 300)
plot_grid(g2, g3, g4, g5, g6, g7, g1)
dev.off()

library(grid)
library(gridExtra)

pdf("F:/Plots/ImputeSingle/fig2e_hrvatin_legend.pdf", 8, 4)
grid.draw(gglegend(g0))
dev.off()











