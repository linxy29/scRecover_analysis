setwd("~/Documents/Data/scRecover/E-MTAB-3929")
load("~/Documents/Data/scRecover/E-MTAB-3929/E_MTAB_3929.Rdata")
counts_original <- E_MTAB_3929[rowSums(E_MTAB_3929 != 0) > 0,]

addTrueZero <- function(path, counts_original){
  print(Sys.time())
  print(paste("Processing files of", path))
  
  scImpute_count <- read.csv(file = paste0(path, "/tempFile/scimpute_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  SAVER_count <- read.csv(file = paste0(path, "/tempFile/SAVER_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  ## whether_impute is the matrix got from ZINB indicating which item should be imputed. A matrix contains true and false.
  ## whether_impute_iz: for non zero value in the original matrix, turn it to false (True are those imputed items); 
  ## whether_impute_inz: for non zero value in the original matrix, turn it to true (True are all non zero values, either it has values in the original matrix or is imputed);
  true_zero = counts_original != 0
  scImute_truezero = scImpute_count * true_zero
  SAVER_truezero = SAVER_count * true_zero
  tempFileDir <- paste0(path, "/tempFile/")
  write.csv(scImute_truezero, file = paste0(tempFileDir, "scImute_truezero.csv"))
  write.csv(SAVER_truezero, file = paste0(tempFileDir, "SAVER_truezero.csv"))
}

summaryImputationSub <- function(path, counts_original, results){
  print(Sys.time())
  print(paste("Processing files of", path))
  
  counts_down <- read.csv(file = paste0(path, "/raw_data.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  scImpute_truezero <- read.csv(file = paste0(path, "/tempFile/scImute_truezero.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  SAVER_truezero <- read.csv(file = paste0(path, "/tempFile/SAVER_truezero.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  
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
  results_scImpute_truezero <- analysisImputation(counts_original, counts_down, scImpute_truezero)
  results_SAVER_truezero <- analysisImputation(counts_original, counts_down, SAVER_truezero)
  results$scImpute_truezero <- results_scImpute_truezero[,3]
  results$SAVER_truezero <- results_SAVER_truezero[,3]
  return(results)
}

load("~/Documents/Data/scRecover/E-MTAB-3929/results_E_MTAB_3929_ImputeSingle_v2.Rdata")
#for(percentage in c(0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.35)){
for(percentage in c(0.25, 0.3, 0.35)){         ## add new results
  for(Kcluster in c(2, 5, 10)){
    for(depth in c(2, 5, 10, 20, 50)){
      resultsName <- paste0("E_MTAB_3929_ImputeSingle_p_", percentage, "_K_", Kcluster, "_d_", depth)
      path <- paste0("./", resultsName)
      #addTrueZero(path, counts_original)
      to_be_saved_obj[[resultsName]] <- summaryImputationSub(path, counts_original, to_be_saved_obj[[resultsName]])
      #assign(resultsName, results)
      #to_be_saved_obj[[resultsName]] <- results
    }
  }
}

save(to_be_saved_obj, file = "~/Documents/Data/scRecover/E-MTAB-3929/results_E_MTAB_3929_ImputeSingle_v3.Rdata")

############### add new statistics
# Function to calculate and add metrics to the data frame
add_metrics_to_df <- function(df) {
  TP <- df["True Positive (Recovered dropout)", ]
  FP <- df["False Positive (Wrong imputed true zero)", ]
  FN <- df["False Negative (Remained dropout)", ]
  TN <- df["True Negative (Remained true zero)", ]
  
  # Calculate the statistics
  df["Precision",] <- TP / (TP + FP)
  #df["Recall",] <- TP / (TP + FN)
  df["Accuracy",] <- (TP + TN) / (TP + FP + FN + TN)
  #df["F1_Score",] <- 2 * (df["Precision",] * df["Recall",]) / (df["Precision",] + df["Recall",])
  df["Specificity",] <- TN / (TN + FP)
  df["False_Positive_Rate",] <- FP / (FP + TN)
  #df["NPV",] <- TN / (TN + FN)
  df["FDR",] <- FP / (TP + FP)
  
  # Calculate MCC
  #mcc_numerator <- TP * TN - FP * FN
  #mcc_denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  #df["MCC",] <- mcc_numerator / mcc_denominator
  
  return(df)
}

# Iterate over each data frame in the list and add the metrics
for (name in names(to_be_saved_obj)) {
  to_be_saved_obj[[name]] <- add_metrics_to_df(to_be_saved_obj[[name]])
}

print("New statsitics like Precision, Accuracy, Specificity, False_Positive_Rate, FDR are added.")

#####################


