setwd("D:/depression/other")
library(dplyr)

file_list <- list.files(pattern = "_rpm_rm0_index_NA.csv")
for (file in file_list) {
  file_path <- file.path("D:/depression/other", file)
  output_file <- gsub("_rpm_rm0_index_NA.csv$", "_de_values_matrix.csv", file) 
  data <- read.csv(file_path)

  control_columns <- grep("^control_", names(data), value = TRUE)
  treatment_columns <- grep("^treatment_", names(data), value = TRUE)

  p_values_matrix <- matrix(NA, nrow = length(data$mirna_name), ncol = 3)
  colnames(p_values_matrix) <- c("mirna_name", "p_value", "Q_value")

  for (i in seq_along(data$mirna_name)) {
    mirna <- data$mirna_name[i]
    control_values <- unlist(data[data$mirna_name == mirna, control_columns])
    treatment_values <- unlist(data[data$mirna_name == mirna, treatment_columns])
    
    result <- try(wilcox.test(treatment_values, control_values))

    if (!inherits(result, "try-error")) {
      p_values_matrix[i, 1] <- mirna
      p_values_matrix[i, 2] <- result$p.value
    } else {
      cat("Error for mirna:", mirna, "\n")
    }
  }
  p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

  print(p_values_matrix)
  write.csv(p_values_matrix, file = file.path("D:/depression/other", output_file), row.names = FALSE)
}