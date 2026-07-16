library(reticulate)
use_condaenv("leiden_env", required = TRUE)

library(readxl)
library(dplyr)
library(stats)
library(ggplot2)
library(factoextra)
library(Rtsne)
library(tibble)
library(tidyr)
library(umap)
library(uwot)
library(RANN)
library(igraph)
library(leiden)
library(caret)
library(nnet)
library(ranger)
library(e1071)  # SVM
library(glmnet) # ElasticNet
library(pROC)   # ROC
library(PRROC)  # PR-Curve
library(kernlab)
library(xgboost)
library(knitr)
library(gridExtra)
library(grid)
library(forcats)

# Set working directory
setwd("/gscratch/geabialn/T_ALL/")

# ==========================================
# 1. READ AND PREPROCESS TRAINING DATA
# ==========================================

# Input expression data (Genes in rows, Patients in columns)
data.raw <- read_excel('data_counts_new.xlsx') %>%
  rename(GeneID = "...1")

print(paste("Original dimensions:", nrow(data.raw), "rows,", ncol(data.raw), "columns"))
as.matrix(data.raw)

# Filter by the Pölönen 300 gene panel
my_genes <- read.table("Polonen_299_genes_new.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]
data.filtered <- data.raw %>% filter(GeneID %in% my_genes)

# Check for missing values
anyNA(data.raw)
print(paste("Number of NAs in raw data:", sum(is.na(data.raw))))

# Save Gene IDs, remove the column to convert the dataframe to numeric
gene_ids <- data.raw$GeneID
expression_matrix <- data.filtered %>%
  select(-GeneID) %>%          
  mutate(across(everything(), as.numeric))
  
# Transpose the data (Samples become rows, Genes become columns)
expression_matrix_t <- as.data.frame(t(expression_matrix))
head(expression_matrix_t)

# Remove zero-variance genes
variance_columns <- apply(expression_matrix_t, 2, var, na.rm = TRUE)
expression_matrix_novariance <- expression_matrix_t[, variance_columns > 0]

# Scale the data (Z-score normalization)
expression_matrix_scaled <- as.data.frame(scale(expression_matrix_novariance))
rownames(expression_matrix_scaled) <- rownames(expression_matrix_novariance)
colnames(expression_matrix_scaled) <- colnames(expression_matrix_novariance)
sum(is.na(expression_matrix_scaled))

# Save scaling parameters (means and standard deviations) for future reuse
scale_means <- attr(scale(expression_matrix_novariance), "scaled:center")
scale_sds <- attr(scale(expression_matrix_novariance), "scaled:scale")
saveRDS(scale_means, "scale_means.rds")
saveRDS(scale_sds, "scale_sds.rds")

# Load clinical subtype labels
labels_raw <- read_excel('T_ALL_labels_full_cohort_total_ordered.xlsx', col_names = TRUE)

# Define features (X) and target variable (y)
X <- expression_matrix_scaled
y <- as.factor(labels_raw$Reviewed.subtype)
sample_ids <- rownames(X) 

# ==========================================
# 2. DATA SPLITTING (STRATIFIED)
# ==========================================
set.seed(54321)
#set.seed(12345)

train_indices <- createDataPartition(
  y = y,
  p = 0.8,
  list = FALSE
)

train_set <- data.frame(X[train_indices, ], y = y[train_indices])
test_set <- data.frame(X[-train_indices, ], y = y[-train_indices])
sample_ids_test <- sample_ids[-train_indices] 

# ==========================================
# 3. DEFINE TUNING GRID & CONTROLS
# ==========================================
xgb_grid <- expand.grid(
  nrounds = 300,
  eta = 0.1,
  max_depth = 6,
  gamma = 0.1,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = c(0.6, 0.8)
)

# Define cross-validation control (5-fold CV)
ctrl <- trainControl(
  method = "cv",
  number = 5,
  verboseIter = TRUE,
  allowParallel = TRUE
)

# ==========================================
# 4. TRAIN THE XGBOOST MODEL & EXTRACT IMPORTANCE
# ==========================================
model_xgb <- train(
  y ~ .,
  data = train_set,
  method = "xgbTree",
  trControl = ctrl,
  tuneGrid = xgb_grid
)


# Save the fully trained model
saveRDS(model_xgb, file = "xgb_model.rds")

cat("\n=== CALCULATING GENE IMPORTANCE ===\n")

# Extract feature importance metrics
importance_xgb <- varImp(model_xgb, scale = TRUE)

# Print the top 20 genes directly to your console
print(importance_xgb, top = 20)

# Generate a visual plot of the top 20 genes
plot(importance_xgb, top = 20, main = "Top 20 Driver Genes for T-ALL Subtyping")

# Export all genes ranked by importance to a CSV file
importance_df <- as.data.frame(importance_xgb$importance)
importance_df <- rownames_to_column(importance_df, "Gene") %>% 
  arrange(desc(Overall))
write.csv(importance_df, "xgb_gene_importance.csv", row.names = FALSE)

# ==========================================
# 5. EVALUATION: CONFUSION MATRIX & PLOT
# ==========================================
predictions_xgb <- predict(model_xgb, newdata = test_set)
cm <- confusionMatrix(predictions_xgb, test_set$y)
cm_table <- as.data.frame(cm$table)

# ------------------------------------------
# PRINT CLEAR METRICS REPORT TO CONSOLE
# ------------------------------------------
cat("\n=======================================================\n")
cat("            MODEL EVALUATION REPORT                    \n")
cat("=======================================================\n")

# 1. Global Metrics
cat(paste(" -> GLOBAL ACCURACY:", round(cm$overall["Accuracy"], 4), "\n"))
cat(paste(" -> 95% Confidence Interval: (", 
          round(cm$overall["AccuracyLower"], 4), ", ", 
          round(cm$overall["AccuracyUpper"], 4), ")\n"))
cat(paste(" -> Cohen's Kappa Index:", round(cm$overall["Kappa"], 4), "\n\n"))

# 2. Extract and Calculate F1-Score per Subtype
metrics_by_class <- as.data.frame(cm$byClass)

# Explicitly calculate F1-Score for each class
metrics_by_class$F1_Score <- (2 * metrics_by_class$Precision * metrics_by_class$Recall) / 
                             (metrics_by_class$Precision + metrics_by_class$Recall)

# Clean up any undefined values (NaN) caused by low prevalence in the test set
metrics_by_class$F1_Score[is.na(metrics_by_class$F1_Score)] <- 0
metrics_by_class$Precision[is.na(metrics_by_class$Precision)] <- 0

# Create a clean summary table with clear metric names
summary_table <- metrics_by_class[, c("Sensitivity", "Specificity", "Precision", "F1_Score")]
colnames(summary_table) <- c("Sensitivity(Recall)", "Specificity", "Precision", "F1-Score")

cat("--- DETAILED METRICS PER SUBTYPE ---\n")
print(round(summary_table, 4))
cat("-------------------------------------------------------\n")

# 3. Average Global Metrics (Macro)
macro_f1 <- mean(summary_table$`F1-Score`, na.rm = TRUE)
macro_precision <- mean(summary_table$Precision, na.rm = TRUE)
macro_recall <- mean(summary_table$`Sensitivity(Recall)`, na.rm = TRUE)

cat(paste(" -> MACRO F1-SCORE AVERAGE :", round(macro_f1, 4), "\n"))
cat(paste(" -> MACRO PRECISION AVERAGE:", round(macro_precision, 4), "\n"))
cat(paste(" -> MACRO RECALL AVERAGE   :", round(macro_recall, 4), "\n"))
cat("=======================================================\n\n")

# Reorder factors for visualization
cm_table <- cm_table %>%
  mutate(
    Reference = fct_reorder(Reference, Freq, .fun = sum, .desc = TRUE),
    Prediction = factor(Prediction, levels = rev(levels(Reference)))
  )

# Plot Confusion Matrix Heatmap
ggplot(cm_table, aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Freq), color = "white") +
  geom_text(aes(label = ifelse(Freq == 0, "", Freq)), size = 4) +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#D1E5F0", "#08306B"),
    values = scales::rescale(c(0, 0.1, max(cm_table$Freq, na.rm = TRUE))),
    na.value = "white",
    name = "Count",
    limits = c(0, max(cm_table$Freq, na.rm = TRUE))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  coord_fixed() +
  labs(
    title = "Confusion Matrix - XGBoost Classifier",
    x = "True Subtype",
    y = "Predicted Subtype"
  )

# ==========================================
# 6. EVALUATION: SAMPLE-WISE PROBABILITIES
# ==========================================
predictions_xgb_prob <- predict(model_xgb, newdata = test_set, type = "prob")
predictions_xgb_prob$TrueSubtype <- test_set$y

threshold <- 0.5

df_prob_long_complete <- predictions_xgb_prob %>%
  pivot_longer(
    cols = -TrueSubtype,
    names_to = "Classifier",
    values_to = "Probability"
  ) %>%
  mutate(
    threshold = 0.5,
    Truth = case_when(
      TrueSubtype == Classifier & Probability >= threshold ~ "True Positive",
      TrueSubtype == Classifier & Probability < threshold ~ "False Negative",
      TrueSubtype != Classifier & Probability >= threshold ~ "False Positive",
      TrueSubtype != Classifier & Probability < threshold ~ "True Negative"
    )
  )

# Plot sample-wise probability distributions
ggplot(df_prob_long_complete, aes(x = Classifier, y = Probability, color = Truth)) +
  geom_jitter(width = 0.25, size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.8, color = "black") +
  scale_color_manual(
    values = c(
      "True Positive" = "#D73027",    # Red
      "False Negative" = "#F46D43",   # Orange
      "False Positive" = "#74ADD1",   # Light Blue
      "True Negative" = "#313695"     # Dark Blue
    )
  ) +
  coord_cartesian(ylim = c(0, 1.02)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right"
  ) +
  labs(
    title = paste0("Sample-wise classifier probabilities (Test cohort, n = ", nrow(test_set), ")"),
    x = "Classifier",
    y = "Predicted probability",
    color = "Prediction Type"
  )

# ==========================================
# 7. IDENTIFYING MISCLASSIFIED SAMPLES
# ==========================================
results_df <- data.frame(
  SampleID = sample_ids_test,
  True_Label = test_set$y,
  Predicted_Label = predictions_xgb,
  stringsAsFactors = FALSE
)

results_df$Correct <- results_df$True_Label == results_df$Predicted_Label

cat("\n=== CLASSIFICATION SUMMARY ===\n")
cat("Total test samples:", nrow(results_df), "\n")
cat("Correctly classified samples:", sum(results_df$Correct), 
    "(", round(mean(results_df$Correct) * 100, 1), "%)\n")
cat("Misclassified samples:", sum(!results_df$Correct), 
    "(", round(mean(!results_df$Correct) * 100, 1), "%)\n\n")

# Audit specific edge case: STAG2/LMO2 misclassified as TME-enriched
stag2_misclassified <- results_df[
  results_df$True_Label == "STAG2/LMO2" & 
  results_df$Predicted_Label == "TME-enriched",
]

# ==========================================
# 8. PIPELINE FOR NEW PATIENT SAMPLE
# ==========================================

# Reload saved model and parameters
model <- readRDS("/beegfs_agamede/gscratch/geabialn/T_ALL/xgb_model.rds")
scale_means <- readRDS("scale_means.rds")
scale_sds <- readRDS("scale_sds.rds")

genes_model_Vnames <- setdiff(colnames(model$trainingData), ".outcome")
genes_original_order <- names(scale_means)

if (length(genes_model_Vnames) != length(genes_original_order)) {
    stop("Error: Number of V-features does not match original gene count.")
}

# Read new sample and handle duplicates
new.raw <- read_excel("counts_new.xlsx") 
names(new.raw)[1] <- "GeneID" 

# Clean duplicates by grouping and summing counts
new.raw_cleaned <- new.raw %>%
  mutate(across(-GeneID, as.numeric)) %>%
  group_by(GeneID) %>%
  summarise(across(where(is.numeric), sum), .groups = 'drop') %>%
  filter(!is.na(GeneID) & GeneID != "") %>%
  mutate(GeneID = as.character(GeneID))

new.raw_matrix <- new.raw_cleaned %>%
  column_to_rownames(var = "GeneID") %>%
  as.matrix()

cat("Dimensions after handling duplicates:", dim(new.raw_matrix), "\n")
cat("Unique gene count:", nrow(new.raw_matrix), "\n")

# Ensure correct sample orientation (genes in rows, sample in column)
if (ncol(new.raw_matrix) > 1) {
  cat("Multiple samples detected. Defaulting to the first column.\n")
  new.raw_matrix <- new.raw_matrix[, 1, drop = FALSE]
}
if (nrow(new.raw_matrix) < ncol(new.raw_matrix)) {
  cat("Transposing matrix to fix orientation (genes to rows)...\n")
  new.raw_matrix <- t(new.raw_matrix)
}

# Clean Ensembl version suffixes
rownames(new.raw_matrix) <- gsub("\\..*", "", rownames(new.raw_matrix))

# Align with the Pölönen gene panel and impute missing ones with 0
genes_present <- intersect(my_genes, rownames(new.raw_matrix))
genes_missing <- setdiff(my_genes, rownames(new.raw_matrix))

cat("Pölönen genes present:", length(genes_present), "\n")
cat("Pölönen genes missing (imputed with 0):", length(genes_missing), "\n")

expr_sub <- new.raw_matrix[genes_present, , drop = FALSE]

if (length(genes_missing) > 0) {
  fill_mat <- matrix(0, 
                     nrow = length(genes_missing), 
                     ncol = ncol(expr_sub),
                     dimnames = list(genes_missing, colnames(expr_sub)))
  expr_sub <- rbind(expr_sub, fill_mat)
}

expr_aligned <- expr_sub[my_genes, , drop = FALSE]

if(any(is.na(expr_aligned))) {
  cat("WARNING: NAs detected. Replacing with 0.\n")
  expr_aligned[is.na(expr_aligned)] <- 0
}

# Rename features to V1...V300 for model structure compatibility
rownames(expr_aligned) <- paste0("V", seq_len(nrow(expr_aligned)))
expression_matrix <- as.data.frame(t(expr_aligned))

# Scale new data using training parameters
cat("Scaling data...\n")
new_expr_scaled <- sweep(expression_matrix, 2, scale_means, "-")
new_expr_scaled <- sweep(new_expr_scaled, 2, scale_sds, "/")
colnames(new_expr_scaled) <- genes_model_Vnames

if(any(is.na(new_expr_scaled))) {
  cat("WARNING: NAs detected after scaling. Fixing...\n")
  na_cols <- colnames(new_expr_scaled)[colSums(is.na(new_expr_scaled)) > 0]
  cat("Columns with NA:", paste(na_cols, collapse=", "), "\n")
  new_expr_scaled[is.na(new_expr_scaled)] <- 0
}

cat("Final dimensions for prediction:", dim(new_expr_scaled), "\n")

# Run prediction
cat("Running prediction...\n")
pred_prob <- predict(model, newdata = new_expr_scaled, type = "prob")
pred_class <- predict(model, newdata = new_expr_scaled)

cat("\n=== PREDICTION RESULTS ===\n")
cat("Predicted Class:", as.character(pred_class), "\n\n")
cat("Probability distribution per class:\n")
print(pred_prob)

# Plot probabilities for the new sample
prob_df <- data.frame(
  Subtype = colnames(pred_prob),
  Probability = as.numeric(pred_prob[1, ])
) %>%
  arrange(desc(Probability))

prob_plot <- ggplot(prob_df, aes(x = reorder(Subtype, Probability), y = Probability)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Prediction Probabilities for New Sample", x = "Subtype", y = "Probability") +
  geom_text(aes(label = sprintf("%.3f", Probability)), hjust = -0.1, size = 3) +
  ylim(0, 1)

print(prob_plot)

# ==========================================
# 9. T-SNE VISUALIZATION WITH NEW SAMPLE
# ==========================================
X_combined <- rbind(expression_matrix_scaled, new_expr_scaled)
y_combined <- c(as.character(labels_raw$Reviewed.subtype), rep("New", nrow(new_expr_scaled)))

set.seed(42)
tsne_result <- Rtsne(X_combined, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 1000)

tsne_df <- data.frame(
  Dim1 = tsne_result$Y[,1],
  Dim2 = tsne_result$Y[,2],
  Subtype = y_combined
)

# Map clinical color palette (adding bright red for the new sample)
colors_subtypes <- c(
  "BCL11B" = "#F8766D", "ETP-like" = "#E7851E", "HOXA9 TCR" = "#D09400",
  "KMT2A" = "#B2A100", "LMO2 gd-like" = "#89AC00", "MLLT10" = "#45B500",
  "NKX2-1" = "#00BC51", "NKX2-5" = "#00C087", "NUP214" = "#00C0B2",
  "NUP98" = "#00C0B3", "SPI1" = "#00BCD6", "STAG2/LMO2" = "#00B3F2",
  "TAL1 aß-like" = "#28A3FB", "TAL1 DP-like" = "#28A3FF", "TLX1" = "#9C8DFF",
  "TLX3" = "#D277FF", "TME-enriched" = "#F166E8", "Others" = "#FF61C7",
  "New" = "red"
)

ggplot(tsne_df, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = Subtype, shape = Subtype), size = 2, alpha = 0.8) +
  scale_color_manual(values = colors_subtypes) +
  scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)) +
  theme_minimal() +
  labs(title = "t-SNE Visualization with New Sample Added")

# ==========================================
# 10. GENERATE AND EXPORT PREDICTION REPORT
# ==========================================
df <- as.data.frame(t(pred_prob)) 
df <- rownames_to_column(df, "Subtype")
df[, -1] <- round(df[, -1], 4)
colnames(df)[2] <- "Prediction"

# Create formatted table object
tabla_grob <- tableGrob(df, rows = NULL, theme = ttheme_default(
  core = list(fg_params = list(cex = 0.6)), 
  colhead = list(fg_params = list(cex = 0.7, fontface = "bold"))
))

# Save report to PDF
ggsave("T_ALL_classifier.pdf", tabla_grob, width = 8, height = 10)