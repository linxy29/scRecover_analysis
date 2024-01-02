library(scRecover)

# Prepare data
#load("F:/data/E-MTAB-3929/E_MTAB_3929.Rdata")
SDRF_original <- read.table("/Users/linxy29/Documents/Data/scRecover/E-MTAB-3929/E-MTAB-3929.sdrf.txt", sep = "\t", header = T, stringsAsFactors = F)
SDRF <- SDRF_original[,c(1, 9, 10, 11, 12, 13, 14, 15)]
row.names(SDRF) <- SDRF[,1]
colnames(SDRF) <- c("sample.name", "individual", "developmental.stage", "treatment", "phenotype", "inferred.lineage", "inferred.TE.subpopulation", "inferred.pseudo.time")

SDRF_E3 <- SDRF[which(SDRF[,"developmental.stage"] == "embryonic day 3"),]
SDRF_E4 <- SDRF[which(SDRF[,"developmental.stage"] == "embryonic day 4"),]
counts_E3E4 <- E_MTAB_3929[, c(SDRF_E3[,"sample.name"], SDRF_E4[,"sample.name"])]
counts_E3E4 <- counts_E3E4[rowSums(counts_E3E4 != 0) > 0,]

# Counts sampling function
countsSampling <- function(counts, fraction){
  n <- floor(sum(counts) * fraction)
  readsGet <- sort(sample(1:sum(counts), n))
  cumCounts <-  c(0, cumsum(counts))
  counts_New <- hist(readsGet, breaks = cumCounts, plot=FALSE)$count
  return(counts_New)
}

# Downsampling 10% reads
set.seed(999)
counts_E3E4_10p <- counts_E3E4
for(i in 1:ncol(counts_E3E4))
  counts_E3E4_10p[,i] <- countsSampling(counts_E3E4[,i], 0.1)
groundTruth <- counts_E3E4 != 0 & counts_E3E4_10p == 0

# Run ImputeSingle
ImputeSingle(counts = counts_E3E4_10p, Kcluster = 2, outputDir = "./E3E4_10p_ImputeSingle/", depth = 20, SAVER = TRUE, MAGIC = TRUE, parallel = TRUE)
scImpute_E3E4 <- read.csv(file = "D:/Research/R_workspace/E3E4_10p_ImputeSingle/tempFile/scimpute_count.csv", stringsAsFactors = FALSE, header = T, row.names = 1)
SAVER_E3E4 <- read.csv(file = "D:/Research/R_workspace/E3E4_10p_ImputeSingle/tempFile/SAVER_count.csv", stringsAsFactors = FALSE, header = T, row.names = 1)
MAGIC_E3E4 <- read.csv(file = "D:/Research/R_workspace/E3E4_10p_ImputeSingle/tempFile/MAGIC_count.csv", stringsAsFactors = FALSE, header = T, row.names = 1)

ImputeSingle_scImpute_E3E4 <- read.csv(file = "D:/Research/R_workspace/E3E4_10p_ImputeSingle/counts_scImpute_iz.csv", stringsAsFactors = FALSE, header = T, row.names = 1)
ImputeSingle_SAVER_E3E4 <- read.csv(file = "D:/Research/R_workspace/E3E4_10p_ImputeSingle/counts_SAVER_iz.csv", stringsAsFactors = FALSE, header = T, row.names = 1)
ImputeSingle_MAGIC_E3E4 <- read.csv(file = "D:/Research/R_workspace/E3E4_10p_ImputeSingle/counts_MAGIC_iz.csv", stringsAsFactors = FALSE, header = T, row.names = 1)

load(file = "./E3E4_10p_ImputeSingle/tempFile/IntermediateVariables.Rdata")

scImputeRight <- groundTruth & (scImpute_E3E4 != 0 & counts_E3E4_10p == 0)
SAVERRight <- groundTruth & (SAVER_E3E4 != 0 & counts_E3E4_10p == 0)
MAGICRight <- groundTruth & (MAGIC_E3E4 != 0 & counts_E3E4_10p == 0)

ImputeSingleRight_scImpute <- groundTruth & (ImputeSingle_scImpute_E3E4 != 0 & counts_E3E4_10p == 0)
ImputeSingleRight_SAVER <- groundTruth & (ImputeSingle_SAVER_E3E4 != 0 & counts_E3E4_10p == 0)
ImputeSingleRight_MAGIC <- groundTruth & (ImputeSingle_MAGIC_E3E4 != 0 & counts_E3E4_10p == 0)

whether_imputeRight <- groundTruth & whether_impute_iz

sum(scImputeRight)
sum(SAVERRight)
sum(MAGICRight)
sum(ImputeSingleRight_scImpute)
sum(ImputeSingleRight_SAVER)
sum(ImputeSingleRight_MAGIC)
sum(whether_imputeRight)

sum(scImpute_E3E4 != 0 & counts_E3E4_10p == 0)
sum(SAVER_E3E4 != 0 & counts_E3E4_10p == 0)
sum(MAGIC_E3E4 != 0 & counts_E3E4_10p == 0)
sum(ImputeSingle_scImpute_E3E4 != 0 & counts_E3E4_10p == 0)
sum(ImputeSingle_SAVER_E3E4 != 0 & counts_E3E4_10p == 0)
sum(ImputeSingle_MAGIC_E3E4 != 0 & counts_E3E4_10p == 0)
sum(whether_impute_iz == TRUE)

sum(scImputeRight) / sum(scImpute_E3E4 != 0 & counts_E3E4_10p == 0)
sum(SAVERRight) / sum(SAVER_E3E4 != 0 & counts_E3E4_10p == 0)
sum(MAGICRight) / sum(MAGIC_E3E4 != 0 & counts_E3E4_10p == 0)
sum(ImputeSingleRight_scImpute) / sum(ImputeSingle_scImpute_E3E4 != 0 & counts_E3E4_10p == 0)
sum(ImputeSingleRight_SAVER) / sum(ImputeSingle_SAVER_E3E4 != 0 & counts_E3E4_10p == 0)
sum(ImputeSingleRight_MAGIC) / sum(ImputeSingle_MAGIC_E3E4 != 0 & counts_E3E4_10p == 0)
sum(whether_imputeRight) / sum(whether_impute_iz == TRUE)
sum(groundTruth) / sum(counts_E3E4_10p == 0)



# Plot imputation results
# PLot Accuracy
load("results_E_MTAB_3929_ImputeSingle.Rdata")
pdf("F:/Plots/ImputeSingle/E_MTAB_3929/Accuracy.pdf", width = 12, height = 9)
par(mfrow=c(3,4))
for(Kcluster in c(2, 5, 10)){
  for(percentage in c(0.01, 0.05, 0.1, 0.2)){
    accuracy_whether_impute <- NULL
    accuracy_scImpute_original <- NULL
    accuracy_scImpute_filter <- NULL
    accuracy_SAVER_original <- NULL
    accuracy_SAVER_filter <- NULL
    accuracy_MAGIC_original <- NULL
    accuracy_MAGIC_filter <- NULL
    for(depth in c(2, 5, 10, 20, 50)){
      resultsName <- paste0("E_MTAB_3929_ImputeSingle_p_", percentage, "_K_", Kcluster, "_d_", depth)
      results <- eval(parse(text = resultsName))
      accuracy_whether_impute <- c(accuracy_whether_impute, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "whether_impute"]*100)
      accuracy_scImpute_original <- c(accuracy_scImpute_original, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "scImpute_original"]*100)
      accuracy_scImpute_filter <- c(accuracy_scImpute_filter, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "scImpute_filter"]*100)
      accuracy_SAVER_original <- c(accuracy_SAVER_original, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "SAVER_original"]*100)
      accuracy_SAVER_filter <- c(accuracy_SAVER_filter, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "SAVER_filter"]*100)
      accuracy_MAGIC_original <- c(accuracy_MAGIC_original, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "MAGIC_original"]*100)
      accuracy_MAGIC_filter <- c(accuracy_MAGIC_filter, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "MAGIC_filter"]*100)
    }
    
    depth <- c(2, 5, 10, 20, 50)
    cex <- 2.5
    lwd <- 2
    allNum <- c(accuracy_SAVER_original, accuracy_scImpute_original, accuracy_SAVER_filter, accuracy_scImpute_filter)
    plot(accuracy_SAVER_original~depth, col="darkolivegreen3", pch=15, lwd=lwd, cex = cex, ylim=c(0, 100), type="b", bty="l", xlab="Depth", ylab="Accuracy / %", xaxt="n", cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.2)
    axis(1, at=depth, labels=depth, cex.axis = 1.2)
    lines(accuracy_MAGIC_original~depth, col="black", pch=18, lwd=lwd, type="b", cex = cex-0.5)
    lines(accuracy_scImpute_original~depth, col="royalblue", pch=17, lwd=lwd, type="b", cex = cex)
    lines(accuracy_SAVER_filter~depth, col="darkolivegreen3", pch=22, lwd=lwd, type="b", cex = cex)
    lines(accuracy_MAGIC_filter~depth, col="black", pch=23, lwd=lwd, type="b", cex = cex-0.8)
    lines(accuracy_scImpute_filter~depth, col="royalblue", pch=24, lwd=lwd, type="b", cex = cex)
    title(paste0("Kcluster = ", Kcluster, "\nDownsampling Rate = ", percentage), cex.main = 1.5)
  }
}
dev.off()

pdf("F:/Plots/ImputeSingle/E_MTAB_3929/Accuracy_legend.pdf", width = 10, height = 3)
plot.new()
legend("center", 
       legend = c("scImpute", "scRecover + scImpute", "SAVER", "scRecover + SAVER", "MAGIC", "scRecover + MAGIC"), 
       col = c("royalblue", "royalblue", "darkolivegreen3", "darkolivegreen3", "black", "black"), 
       pch = c(17, 24, 15, 22, 18, 23), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.1,
       lwd = 2,
       ncol = 3,
       pt.bg = 'white',
       text.col = "black",
       horiz = F)
dev.off()



# PLot predicted dropout number
load("results_E_MTAB_3929_ImputeSingle.Rdata")
pdf("F:/Plots/ImputeSingle/E_MTAB_3929/Dropout.pdf", width = 12, height = 9)
par(mfrow=c(3,4))
for(Kcluster in c(2, 5, 10)){
  for(percentage in c(0.01, 0.05, 0.1, 0.2)){
    dropout_downsampling <- NULL
    dropout_whether_impute <- NULL
    dropout_scImpute_original <- NULL
    dropout_scImpute_filter <- NULL
    dropout_SAVER_original <- NULL
    dropout_SAVER_filter <- NULL
    dropout_MAGIC_original <- NULL
    dropout_MAGIC_filter <- NULL
    for(depth in c(2, 5, 10, 20, 50)){
      resultsName <- paste0("E_MTAB_3929_ImputeSingle_p_", percentage, "_K_", Kcluster, "_d_", depth)
      results <- eval(parse(text = resultsName))
      dropout_downsampling <- c(dropout_downsampling, results["Dropout number per cell", "downsampling"])
      dropout_whether_impute <- c(dropout_whether_impute, results["Imputed zero to non-zero number", "whether_impute"])
      dropout_scImpute_original <- c(dropout_scImpute_original, results["Imputed zero to non-zero number", "scImpute_original"])
      dropout_scImpute_filter <- c(dropout_scImpute_filter, results["Imputed zero to non-zero number", "scImpute_filter"])
      dropout_SAVER_original <- c(dropout_SAVER_original, results["Imputed zero to non-zero number", "SAVER_original"])
      dropout_SAVER_filter <- c(dropout_SAVER_filter, results["Imputed zero to non-zero number", "SAVER_filter"])
      dropout_MAGIC_original <- c(dropout_MAGIC_original, results["Imputed zero to non-zero number", "MAGIC_original"])
      dropout_MAGIC_filter <- c(dropout_MAGIC_filter, results["Imputed zero to non-zero number", "MAGIC_filter"])
    }
    
    depth <- c(2, 5, 10, 20, 50)
    cex <- 2.5
    lwd <- 2
    allNum <- c(dropout_downsampling, dropout_SAVER_original, dropout_scImpute_original, dropout_SAVER_filter, dropout_scImpute_filter)
    plot(dropout_SAVER_original~depth, col="darkolivegreen3", pch=15, lwd=lwd, cex = cex, ylim=c(0, 16000), type="b", bty="l", xlab="Depth", ylab="Predicted dropout number", xaxt="n", cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.2)
    axis(1, at=depth, labels=depth, cex.axis = 1.2)
    lines(dropout_MAGIC_original~depth, col="black", pch=18, lwd=lwd, type="b", cex = cex-0.5)
    lines(dropout_scImpute_original~depth, col="royalblue", pch=17, lwd=lwd, type="b", cex = cex)
    lines(dropout_SAVER_filter~depth, col="darkolivegreen3", pch=22, lwd=lwd, type="b", cex = cex)
    lines(dropout_MAGIC_filter~depth, col="black", pch=23, lwd=lwd, type="b", cex = cex-0.8)
    lines(dropout_scImpute_filter~depth, col="royalblue", pch=24, lwd=lwd, type="b", cex = cex)
    lines(dropout_downsampling~depth, col="red", pch=20, lwd=lwd, type="b", cex = cex)
    title(paste0("Kcluster = ", Kcluster, "\nDownsampling Rate = ", percentage), cex.main = 1.5)
  }
}
dev.off()

pdf("F:/Plots/ImputeSingle/E_MTAB_3929/Dropout_legend.pdf", width = 12, height = 3)
plot.new()
legend("center", 
       legend = c("scImpute", "scRecover + scImpute", "SAVER", "scRecover + SAVER", "MAGIC", "scRecover + MAGIC", "Real dropout"), 
       col = c("royalblue", "royalblue", "darkolivegreen3", "darkolivegreen3", "black", "black", "red"), 
       pch = c(17, 24, 15, 22, 18, 23, 20), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.1,
       lwd = 2,
       ncol = 4,
       pt.bg = 'white',
       text.col = "black",
       horiz = F)
dev.off()










# Plot imputation results KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK = 1
# PLot Accuracy and predicted dropout number
# 中文显示在plot后面加了family = 'SimSun'，如果不用中文则须去掉
library(showtext)
showtext.auto(enable = TRUE)
font.add('SimSun', 'simsun.ttc')
# 中文显示在plot后面加了family = 'SimSun'，如果不用中文则须去掉

load("results_E_MTAB_3929_ImputeSingle_K_1.Rdata")
pdf("F:/Plots/ImputeSingle/E_MTAB_3929/K1_Accuracy_Dropout.pdf", width = 12, height = 6)
# png("F:/Plots/ImputeSingle/E_MTAB_3929/K1_Accuracy_Dropout.png", width = 1000, height = 500)
cex <- 2.5
lwd <- 2
par(mfrow=c(2,4))
for(Kcluster in c(1)){
  for(percentage in c(0.01, 0.05, 0.1, 0.2)){
    accuracy_whether_impute <- NULL
    accuracy_scImpute_original <- NULL
    accuracy_scImpute_filter <- NULL
    accuracy_SAVER_original <- NULL
    accuracy_SAVER_filter <- NULL
    accuracy_MAGIC_original <- NULL
    accuracy_MAGIC_filter <- NULL
    for(depth in c(2, 5, 10, 20, 50)){
      resultsName <- paste0("E_MTAB_3929_ImputeSingle_p_", percentage, "_K_", Kcluster, "_d_", depth)
      results <- eval(parse(text = resultsName))
      accuracy_whether_impute <- c(accuracy_whether_impute, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "whether_impute"]*100)
      accuracy_scImpute_original <- c(accuracy_scImpute_original, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "scImpute_original"]*100)
      accuracy_scImpute_filter <- c(accuracy_scImpute_filter, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "scImpute_filter"]*100)
      accuracy_SAVER_original <- c(accuracy_SAVER_original, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "SAVER_original"]*100)
      accuracy_SAVER_filter <- c(accuracy_SAVER_filter, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "SAVER_filter"]*100)
      accuracy_MAGIC_original <- c(accuracy_MAGIC_original, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "MAGIC_original"]*100)
      accuracy_MAGIC_filter <- c(accuracy_MAGIC_filter, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "MAGIC_filter"]*100)
    }
    
    depth <- c(2, 5, 10, 20, 50)
    allNum <- c(accuracy_SAVER_original, accuracy_scImpute_original, accuracy_SAVER_filter, accuracy_scImpute_filter)
    plot(accuracy_SAVER_original~depth, col="darkolivegreen3", pch=15, lwd=lwd, cex = cex, ylim=c(0, 100), type="b", bty="l", xlab="预测深度", ylab="准确率 / %", xaxt="n", cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.2, family = 'SimSun')
    axis(1, at=depth, labels=depth, cex.axis = 1.2)
    lines(accuracy_MAGIC_original~depth, col="black", pch=18, lwd=lwd, type="b", cex = cex-0.5)
    lines(accuracy_scImpute_original~depth, col="royalblue", pch=17, lwd=lwd, type="b", cex = cex)
    lines(accuracy_SAVER_filter~depth, col="darkolivegreen3", pch=22, lwd=lwd, type="b", cex = cex)
    lines(accuracy_MAGIC_filter~depth, col="black", pch=23, lwd=lwd, type="b", cex = cex-0.8)
    lines(accuracy_scImpute_filter~depth, col="royalblue", pch=24, lwd=lwd, type="b", cex = cex)
    title(paste0("降采样率 = ", percentage), cex.main = 1.5)
  }
}

for(Kcluster in c(1)){
  for(percentage in c(0.01, 0.05, 0.1, 0.2)){
    dropout_downsampling <- NULL
    dropout_whether_impute <- NULL
    dropout_scImpute_original <- NULL
    dropout_scImpute_filter <- NULL
    dropout_SAVER_original <- NULL
    dropout_SAVER_filter <- NULL
    dropout_MAGIC_original <- NULL
    dropout_MAGIC_filter <- NULL
    for(depth in c(2, 5, 10, 20, 50)){
      resultsName <- paste0("E_MTAB_3929_ImputeSingle_p_", percentage, "_K_", Kcluster, "_d_", depth)
      results <- eval(parse(text = resultsName))
      dropout_downsampling <- c(dropout_downsampling, results["Dropout number per cell", "downsampling"])
      dropout_whether_impute <- c(dropout_whether_impute, results["Imputed zero to non-zero number", "whether_impute"])
      dropout_scImpute_original <- c(dropout_scImpute_original, results["Imputed zero to non-zero number", "scImpute_original"])
      dropout_scImpute_filter <- c(dropout_scImpute_filter, results["Imputed zero to non-zero number", "scImpute_filter"])
      dropout_SAVER_original <- c(dropout_SAVER_original, results["Imputed zero to non-zero number", "SAVER_original"])
      dropout_SAVER_filter <- c(dropout_SAVER_filter, results["Imputed zero to non-zero number", "SAVER_filter"])
      dropout_MAGIC_original <- c(dropout_MAGIC_original, results["Imputed zero to non-zero number", "MAGIC_original"])
      dropout_MAGIC_filter <- c(dropout_MAGIC_filter, results["Imputed zero to non-zero number", "MAGIC_filter"])
    }
    
    depth <- c(2, 5, 10, 20, 50)
    allNum <- c(dropout_downsampling, dropout_SAVER_original, dropout_scImpute_original, dropout_SAVER_filter, dropout_scImpute_filter)
    plot(dropout_SAVER_original~depth, col="darkolivegreen3", pch=15, lwd=lwd, cex = cex, ylim=c(0, 16000), type="b", bty="l", xlab="预测深度", ylab="填补漏失零数", xaxt="n", cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.2, family = 'SimSun')
    axis(1, at=depth, labels=depth, cex.axis = 1.2)
    lines(dropout_MAGIC_original~depth, col="black", pch=18, lwd=lwd, type="b", cex = cex-0.5)
    lines(dropout_scImpute_original~depth, col="royalblue", pch=17, lwd=lwd, type="b", cex = cex)
    lines(dropout_SAVER_filter~depth, col="darkolivegreen3", pch=22, lwd=lwd, type="b", cex = cex)
    lines(dropout_MAGIC_filter~depth, col="black", pch=23, lwd=lwd, type="b", cex = cex-0.8)
    lines(dropout_scImpute_filter~depth, col="royalblue", pch=24, lwd=lwd, type="b", cex = cex)
    lines(dropout_downsampling~depth, col="red", pch=20, lwd=lwd, type="b", cex = cex)
    title(paste0("降采样率 = ", percentage), cex.main = 1.5)
  }
}
dev.off()

pdf("F:/Plots/ImputeSingle/E_MTAB_3929/K1_legend.pdf", width = 12, height = 3)
plot.new()
legend("center", 
       legend = c("scImpute", "scRecover + scImpute", "SAVER", "scRecover + SAVER", "MAGIC", "scRecover + MAGIC", "Real dropout"), 
       col = c("royalblue", "royalblue", "darkolivegreen3", "darkolivegreen3", "black", "black", "red"), 
       pch = c(17, 24, 15, 22, 18, 23, 20), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.1,
       lwd = 2,
       ncol = 4,
       pt.bg = 'white',
       text.col = "black",
       horiz = F)
dev.off()





# Overlap of scRecover imputation
path <- "F:/Plots/ImputeSingle/E_MTAB_3929/E_MTAB_3929_ImputeSingle_p_0.1_K_1_d_20/"
K1 <- readRDS(file = paste0(path, "whether_impute_iz_p_0.1_K_1_d_20.rds"))
K2 <- readRDS(file = paste0(path, "whether_impute_iz_p_0.1_K_2_d_20.rds"))
K5 <- readRDS(file = paste0(path, "whether_impute_iz_p_0.1_K_5_d_20.rds"))
K10 <- readRDS(file = paste0(path, "whether_impute_iz_p_0.1_K_10_d_20.rds"))
whether_impute_iz <- list(K1 = as.matrix(K1), K2 = as.matrix(K2), K5 = as.matrix(K5), K10 = as.matrix(K10))
overlap <- matrix(NA, nrow = 4, ncol = 4, dimnames = list(names(whether_impute_iz), names(whether_impute_iz)))
for(i in 1:4){
  for(j in 1:4){
    overlap[i, j] <- round(mean(colSums(whether_impute_iz[[i]] & whether_impute_iz[[j]])))
  }
}





# Heatmap of downsampling and imputed data
library(gplots)
library(scales)
load("F:/data/E-MTAB-3929/E_MTAB_3929.Rdata")
path <- "F:/Plots/ImputeSingle/E_MTAB_3929/E_MTAB_3929_ImputeSingle_p_0.1_K_1_d_20/"
counts_original <- E_MTAB_3929[rowSums(E_MTAB_3929 != 0) > 0,]
counts_down <- readRDS(file = paste0(path, "raw_data.rds"))
scImpute_filter <- readRDS(file = paste0(path, "scImpute_filter.rds"))
SAVER_filter <- readRDS(file = paste0(path, "SAVER_filter.rds"))
MAGIC_filter <- readRDS(file = paste0(path, "MAGIC_filter.rds"))
scimpute_count <- readRDS(file = paste0(path, "tempFile/scimpute_count.rds"))
SAVER_count <- readRDS(file = paste0(path, "tempFile/SAVER_count.rds"))
MAGIC_count <- readRDS(file = paste0(path, "tempFile/MAGIC_count.rds"))


sampleNum <- ncol(counts_original)
counts_family_SAVER <- as.matrix(cbind(counts_original, counts_down, SAVER_count, SAVER_filter))
counts_family_SAVER <- counts_family_SAVER[rowSums(counts_down != 0) > 0,]
counts_family_SAVER <- log(counts_family_SAVER + 1)
counts_family_SAVER <- counts_family_SAVER[sample(1:nrow(counts_family_SAVER), 1000),]

mainText <- "counts family of SAVER"
pdf(file = "D:/temp.pdf", width = 8, height = 8)
heatmap.2(counts_family_SAVER, col = c("blue", rev(heat.colors(1024))), lhei = c(1,5.7), ColSideColors = rep(c(hue_pal()(4)[1], hue_pal()(4)[2], hue_pal()(4)[3], hue_pal()(4)[4]), each = sampleNum), dendrogram = "row", Rowv = T, Colv = F, scale = "none", labRow = NA, labCol = NA, trace="none", main = mainText)
legend(0.32, 1.04, legend = c("Original", "Downsampling", "SAVER", "scRecover + SAVER"), horiz = T, col = c(hue_pal()(4)[1], hue_pal()(4)[2], hue_pal()(4)[3], hue_pal()(4)[4]), pch = 15, pt.cex=2, xpd=TRUE)
dev.off()









