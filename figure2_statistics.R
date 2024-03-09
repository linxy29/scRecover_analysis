#### figure 2a

# Install and load necessary libraries
library(tidyverse)

setwd("~/Documents/Data/scRecover/E-MTAB-3929")

# load result file
load("~/Documents/Data/scRecover/E-MTAB-3929/results_E_MTAB_3929_ImputeSingle_v2.Rdata")

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

############### Plot the results (supplementary)
# Define the metrics names to extract
metrics <- c("Precision", "Accuracy", "Specificity", "False_Positive_Rate", "FDR")

# Create a list to hold each metric's data frames
metrics_data <- setNames(vector("list", length(metrics)), metrics)

# Iterate through each metric and extract the corresponding data across all data frames
for (metric in metrics) {
  # Extract the metric data and the name of the data frame and store in a list
  metrics_data[[metric]] <- mapply(function(df, name) {
    if (metric %in% rownames(df)) {
      return(data.frame(Name = name, Value = df[metric, ], stringsAsFactors = FALSE))
    } else {
      return(data.frame(Name = name, Value = NA, stringsAsFactors = FALSE)) # In case the metric is not found
    }
  }, to_be_saved_obj, names(to_be_saved_obj), SIMPLIFY = FALSE)
  
  # Combine the list of data frames for the current metric into one data frame
  metrics_data[[metric]] <- do.call(rbind, metrics_data[[metric]])
}

# At this point, 'metrics_data' is a list where each element is a data frame for a specific metric
# Each data frame contains the metric values and the names of the original data frames

# plot 
library(ggplot2)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# Assuming each data frame in metrics_data list has the same structure as your example for Precision

# Define a function to process each data frame
process_and_plot <- function(df, metric_name) {
  # 1. Extract p, K, d from the Name column
  data <- df %>%
    mutate(
      p = as.numeric(str_extract(Name, "(?<=_p_)[0-9.]+")),
      K = as.numeric(str_extract(Name, "(?<=_K_)[0-9]+")),
      d = as.numeric(str_extract(Name, "(?<=_d_)[0-9]+"))) %>%
    select(-Name, -Value.original, -Value.downsampling) %>%
    pivot_longer(
      cols = starts_with("Value"),
      names_to = "Method",
      values_to = "Value"
    ) %>%
    drop_na(Value)
  
  # 2. Group by p and d, and select the best K for each metric
  # Determine the function to use based on the metric type
  best_value_function <- case_when(
    metric_name %in% c("Precision", "Accuracy", "Specificity") ~ "max",
    metric_name %in% c("False_Positive_Rate", "FDR") ~ "min",
    TRUE ~ NA_character_  # This should not happen but provides a fallback
  )
  
  # Add a column to indicate the best value row
  data <- data %>%
    group_by(p, d, Method) %>%
    mutate(
      Best_Value = case_when(
        best_value_function == "max" & Value == max(Value, na.rm = TRUE) ~ TRUE,
        best_value_function == "min" & Value == min(Value, na.rm = TRUE) ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    filter(Best_Value) %>%
    filter(d >= 10) %>% 
    mutate(Method = str_remove(Method, "Value.")) %>%
    filter(Method %in% c("SAVER_filter", "scImpute_filter", "SAVER_original", "scImpute_original")) %>%
    mutate(Method = str_replace(Method, "_filter", " + scRecover")) %>% 
    filter(p != 0.01) %>% 
    ungroup()
  
  # 3. Plot the change of values for different method when the d change
  plot <- ggplot(data, aes(x = d, y = Value, color = Method)) +
    geom_line() +
    geom_point() +
    facet_wrap(~p, scales = "free") +
    labs(title = paste("Change of", metric_name, "for different methods and d values"),
         x = "Prediction Depth (d)",
         y = metric_name) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 4. Save the plot to a file
  # Create a filename based on the metric name
  file_name <- paste0("summary_plot/plot_", metric_name, ".pdf")
  ggsave(file_name, plot, width = 10, height = 6, dpi = 300)
  
  return(data)
}

# Initialize an empty list to store processed data frames
processed_dfs <- list()

# Iterate over each metric, process the data, and store the result with metric_name
for(metric_name in names(metrics_data)) {
  print(paste("Processing", metric_name))
  df <- metrics_data[[metric_name]]
  processed_dfs[[metric_name]] <- process_and_plot(df, metric_name)
}

# Save the list of processed data frames 
saveRDS(processed_dfs, "summary_statistics_frames.rds")

# Load the openxlsx package
library(openxlsx)
# Create a new workbook
wb <- createWorkbook()
# Add sheets for each item in processed_dfs list
for (name in names(processed_dfs)) {
  # Create a new sheet with the name of the dataframe
  addWorksheet(wb, name)
  
  # Write the dataframe to the newly created sheet
  writeData(wb, name, processed_dfs[[name]])
}
# Save the workbook to a file
saveWorkbook(wb, "summary_statistics_frames.xlsx", overwrite = TRUE)
# Print message to indicate completion
cat("Excel file 'summary_statistics_frames.xlsx' has been created with separate sheets for Precision, Recall, and Accuracy.")



############### Plot the results figure 2

metrics_of_interest <- c("Accuracy", "Precision", "Specificity")
p_values_of_interest <- c(0.05, 0.15, 0.25, 0.35)

# Combine data frames into one, preserving metric names
combined_df <- bind_rows(processed_dfs[metrics_of_interest], .id = "Metric")

# Adjusting facet_wrap to organize subplots by both Metric and p, in different columns
combined_df %>% 
  filter(p %in% p_values_of_interest) %>% 
  mutate(p = str_c("p = ", p)) %>% 
ggplot(aes(x = d, y = Value, color = Method)) +
  geom_line() +
  geom_point() +
  facet_grid(Metric ~ p, scales = "free_y") +  # Changed to facet_grid for better control
  labs(title = "Metrics Comparison Across Different Scenarios",
       x = "Prediction Depth (d)",
       y = "Metric Value") +
  theme_classic() +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 0))  # Adjust text angle if needed

ggsave("figure2.pdf")
