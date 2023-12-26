library(ImputeSingle)
setwd("/data/miaozhun/analysis/imputation/Chu1")
counts_original <- readRDS(file = "/data/miaozhun/analysis/imputation/Chu1/Chu1.rds")
counts_original <- as.matrix(counts_original)

# Counts sampling function
countsSampling <- function(counts, fraction){
  n <- floor(sum(counts) * fraction)
  readsGet <- sort(sample(1:sum(counts), n))
  cumCounts <-  c(0, cumsum(counts))
  counts_New <- hist(readsGet, breaks = cumCounts, plot=FALSE)$count
  return(counts_New)
}

# Get dropout matrix function
getDropoutMatrix <- function(counts){
  lambda <- 0.1
  counts <- as.matrix(counts)
  dropoutMatrix <- counts
  dropoutMatrix[,] <- 0
  dropoutRateAll <- NULL
  for(i in 1:nrow(counts)){
    cat("\r",paste0("Generating ", i, " of ", nrow(counts), " genes with dropouts..."))
    nonzero <- counts[i,][counts[i,] != 0]
    dropoutRate <- exp(-lambda * (mean(log(nonzero))^2))
    if(is.nan(dropoutRate))
      dropoutRate <- 1
    dropoutRateAll <- c(dropoutRateAll, dropoutRate)
    dropoutMatrix[i,] <- rbinom(ncol(counts), 1, 1 - dropoutRate)
  }
  write(paste0("Mean dropoutRate = ", round(mean(dropoutRateAll), digits = 3)), file = "./dropoutRate.txt", append = TRUE)
  write(quantile(dropoutRateAll), file = "./dropoutRate.txt", append = TRUE)
  write(paste0("Mean dropoutRate(!= 1) = ", round(mean(dropoutRateAll[dropoutRateAll != 1]), digits = 3)), file = "./dropoutRate.txt", append = TRUE)
  write(quantile(dropoutRateAll[dropoutRateAll != 1]), file = "./dropoutRate.txt", append = TRUE)
  write("\n", file = "./dropoutRate.txt", append = TRUE)
  return(dropoutMatrix)
}

for(percentage in c(0.01, 0.05, 0.1, 0.2)){
  # Downsampling reads
  set.seed(999)
  counts_down <- counts_original
  for(i in 1:ncol(counts_original))
    counts_down[,i] <- countsSampling(counts_original[,i], percentage)
  dropoutMatrix <- getDropoutMatrix(counts_original)
  counts_down <- counts_down * dropoutMatrix
  
  # Run ImputeSingle
  for(Kcluster in c(2, 5, 10)){
    for(depth in c(2, 5, 10, 20, 50)){
      outputDir <- paste0("./Chu1_ImputeSingle_p_", percentage, "_K_", Kcluster, "_d_", depth, "/")
      ImputeSingle(counts = counts_down, Kcluster = Kcluster, outputDir = outputDir, depth = depth, SAVER = TRUE, MAGIC = TRUE, parallel = TRUE)
      gc()
    }
  }
}





# Summary Imputation results
summaryImputation <- function(path, counts_original){
  print(Sys.time())
  print(paste("Processing files of", path))
  
  counts_down <- read.csv(file = paste0(path, "/raw_data.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  whether_impute_inz <- readRDS(file = paste0(path, "/tempFile/whether_impute_inz.rds"))
  scImpute_count <- read.csv(file = paste0(path, "/tempFile/scimpute_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  SAVER_count <- read.csv(file = paste0(path, "/tempFile/SAVER_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  MAGIC_count <- read.csv(file = paste0(path, "/tempFile/MAGIC_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  scImpute_filter <- read.csv(file = paste0(path, "/scImpute_filter.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  SAVER_filter <- read.csv(file = paste0(path, "/SAVER_filter.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  MAGIC_filter <- read.csv(file = paste0(path, "/MAGIC_filter.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  
  analysisImputation <- function(counts_original, counts_down, counts_imputation){
    feature <- c("Expressed gene number", "Non-zero number per cell", "Zero number per cell")
    results <- data.frame(original = rep(NA, times = length(feature)), downsampling = rep(NA, times = length(feature)), imputation = rep(NA, times = length(feature)), row.names = feature)
    results["Expressed gene number", ] <- c(sum(rowSums(counts_original != 0) > 0), sum(rowSums(counts_down != 0) > 0), sum(rowSums(counts_imputation != 0) > 0))
    results["Non-zero number per cell", ] <- c(mean(colSums(counts_original != 0)), mean(colSums(counts_down != 0)), mean(colSums(counts_imputation != 0)))
    results["Zero number per cell", ] <- c(mean(colSums(counts_original == 0)), mean(colSums(counts_down == 0)), mean(colSums(counts_imputation == 0)))
    results["Dropout number per cell",] <- c(0, mean(colSums(counts_down == 0 & counts_original != 0)), mean(colSums(counts_imputation == 0 & counts_original != 0)))
    results["Imputed zero to non-zero number",] <- c(NA, NA, mean(colSums(counts_imputation != 0 & counts_down == 0)))
    results["True Positive (Recovered dropout)",] <- c(NA, NA, mean(colSums(counts_imputation != 0 & counts_down == 0 & counts_original != 0)))
    results["False Positive (Wrong imputed true zero)",] <- c(NA, NA, mean(colSums(counts_imputation != 0 & counts_down == 0 & counts_original == 0)))
    results["False Negative (Remained dropout)",] <- c(NA, NA, mean(colSums(counts_imputation == 0 & counts_original != 0)))
    results["True Negative (Remained true zero)",] <- c(NA, NA, mean(colSums(counts_imputation == 0 & counts_original == 0)))
    results["Precision = TP / (TP + FP)",] <- c(NA, NA, results["True Positive (Recovered dropout)",3]/(results["True Positive (Recovered dropout)",3] + results["False Positive (Wrong imputed true zero)",3]))
    results["Recall = TP / (TP + FN)",] <- c(NA, NA, results["True Positive (Recovered dropout)",3]/(results["True Positive (Recovered dropout)",3] + results["False Negative (Remained dropout)",3]))
    results["Accuracy = (TP + TN) / (TP + FP + TN + FN)",] <- c(NA, NA, (results["True Positive (Recovered dropout)",3] + results["True Negative (Remained true zero)",3])/(results["True Positive (Recovered dropout)",3] + results["False Positive (Wrong imputed true zero)",3] + results["False Negative (Remained dropout)",3] + results["True Negative (Remained true zero)",3]))
    results <- round(results, digits = 3)
    return(results)
  }
  
  results_whether_impute_inz <- analysisImputation(counts_original, counts_down, whether_impute_inz)
  results_scImpute_count <- analysisImputation(counts_original, counts_down, scImpute_count)
  results_scImpute_filter <- analysisImputation(counts_original, counts_down, scImpute_filter)
  results_SAVER_count <- analysisImputation(counts_original, counts_down, SAVER_count)
  results_SAVER_filter <- analysisImputation(counts_original, counts_down, SAVER_filter)
  results_MAGIC_count <- analysisImputation(counts_original, counts_down, MAGIC_count)
  results_MAGIC_filter <- analysisImputation(counts_original, counts_down, MAGIC_filter)
  results <- cbind(results_whether_impute_inz, results_scImpute_count[,3], results_scImpute_filter[,3], results_SAVER_count[,3], results_SAVER_filter[,3], results_MAGIC_count[,3], results_MAGIC_filter[,3])
  colnames(results)[3:9] <- c("whether_impute", "scImpute_original", "scImpute_filter", "SAVER_original", "SAVER_filter", "MAGIC_original", "MAGIC_filter")
  return(results)
}

library(Matrix)
setwd("/data/miaozhun/analysis/imputation/Chu1")
counts_original <- readRDS(file = "/data/miaozhun/analysis/imputation/Chu1/Chu1.rds")
counts_original <- as.matrix(counts_original)
to_be_saved_obj <- NULL

for(percentage in c(0.01, 0.05, 0.1, 0.2)){
  for(Kcluster in c(2, 5, 10)){
    for(depth in c(2, 5, 10, 20, 50)){
      resultsName <- paste0("Chu1_ImputeSingle_p_", percentage, "_K_", Kcluster, "_d_", depth)
      path <- paste0("./", resultsName)
      results <- summaryImputation(path, counts_original)
      assign(resultsName, results)
      to_be_saved_obj <- c(to_be_saved_obj, resultsName)
    }
    save(list = to_be_saved_obj, file = "./results_Chu1_ImputeSingle.Rdata")
  }
}

save(list = to_be_saved_obj, file = "./results_Chu1_ImputeSingle.Rdata")




