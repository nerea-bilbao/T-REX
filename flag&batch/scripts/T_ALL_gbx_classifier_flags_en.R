# ============================================================
# 1. LIBRARIES
# ============================================================

library(readxl)
library(dplyr)
library(caret)
library(e1071)
library(xgboost)
library(pROC)
library(ggplot2)
library(tidyr)
library(tibble)
library(Rtsne)
library(gridExtra)

set.seed(54321)

# ============================================================
# 2. LOADING AND PREPROCESSING
# ============================================================

data.raw <- read_excel('data_counts_new.xlsx') %>%
  rename(GeneID = "...1")

my_genes <- read.table("Polonen_299_genes_new.txt",
                       header = FALSE,
                       stringsAsFactors = FALSE)[[1]]

data.filtered <- data.raw %>% filter(GeneID %in% my_genes)

expression_matrix <- data.filtered %>%
  select(-GeneID) %>%
  mutate(across(everything(), as.numeric))

expression_matrix_t <- as.data.frame(t(expression_matrix))

variance_columns <- apply(expression_matrix_t, 2, var)
expression_matrix_novariance <- expression_matrix_t[, variance_columns > 0]

expression_matrix_scaled <- as.data.frame(scale(expression_matrix_novariance))

scale_means <- attr(scale(expression_matrix_novariance), "scaled:center")
scale_sds   <- attr(scale(expression_matrix_novariance), "scaled:scale")

saveRDS(scale_means, "scale_means.rds")
saveRDS(scale_sds, "scale_sds.rds")

labels_raw <- read_excel('T_ALL_labels_full_cohort_total_ordered.xlsx')

X <- expression_matrix_scaled
y_original <- as.factor(labels_raw$Reviewed.subtype)
y <- y_original
levels(y) <- make.names(levels(y), unique = TRUE)

# Save mapping to restore original names in predictions
name_mapping <- setNames(levels(y_original), levels(y))
saveRDS(name_mapping, "subtype_name_mapping.rds")


# ============================================================
# 3. TRAIN / TEST SPLIT
# ============================================================
set.seed(12345)

train_idx <- createDataPartition(y, p = 0.8, list = FALSE)

X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]
y_train <- y[train_idx]
y_test  <- y[-train_idx]

train_data <- data.frame(X_train, y = y_train)
test_data  <- data.frame(X_test,  y = y_test)

# ============================================================
# 4. XGBOOST TRAINING
# ============================================================

xgb_grid <- expand.grid(
  nrounds = 300,
  eta = 0.1,
  max_depth = 6,
  gamma = 0.1,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)

ctrl <- trainControl(method = "cv",
                     number = 5,
                     classProbs = TRUE)

model_xgb <- train(
  y ~ .,
  data = train_data,
  method = "xgbTree",
  trControl = ctrl,
  tuneGrid = xgb_grid,
  metric = "Kappa"
)

saveRDS(model_xgb, "xgb_model_full.rds")

# ============================================================
# 5. GLOBAL OOD (OUT-OF-DISTRIBUTION) DETECTOR
# ============================================================

ood_detector <- svm(
  X_train,
  type = "one-classification",
  nu = 0.1,
  gamma = 1/ncol(X_train),
  scale = FALSE
)

saveRDS(ood_detector, "ood_detector_full.rds")

# ============================================================
# 6. UNCERTAINTY METRICS ON TEST SET
# ============================================================

test_probs <- predict(model_xgb, newdata = test_data, type = "prob")
test_pred  <- predict(model_xgb, newdata = test_data)

test_max_prob <- apply(test_probs, 1, max)

test_entropy <- apply(test_probs, 1, function(p)
  -sum(p * log(p + 1e-12))
)

test_correct <- test_pred == y_test

uncertainty_df <- data.frame(
  max_prob = test_max_prob,
  entropy  = test_entropy,
  correct  = test_correct
)

# ============================================================
# 7. ROC OPTIMIZATION OF THRESHOLDS
# ============================================================

roc_prob <- roc(test_correct, test_max_prob)
prob_threshold <- as.numeric(coords(roc_prob, "best",
                                    best.method = "youden",
                                    ret = "threshold")$threshold)

roc_entropy <- roc(test_correct, -test_entropy)
entropy_threshold <- as.numeric(coords(roc_entropy, "best",
                                       best.method = "youden",
                                       ret = "threshold")$threshold)


saveRDS(prob_threshold, "prob_threshold_full.rds")
saveRDS(entropy_threshold, "entropy_threshold_full.rds")

cat("Probability threshold:", prob_threshold, "\n")
cat("Entropy threshold:", entropy_threshold, "\n")

# ============================================================
# 8. COVERAGE vs RESIDUAL ERROR EVALUATION
# ============================================================

prob_threshold <- as.numeric(prob_threshold)
entropy_threshold <- as.numeric(entropy_threshold)

flag_vector <- (uncertainty_df$max_prob < prob_threshold) |
  (uncertainty_df$entropy > entropy_threshold)

coverage <- mean(!flag_vector)
residual_error <- mean(!uncertainty_df$correct[!flag_vector])

cat("Coverage:", coverage, "\n")
cat("Residual error after filtering:", residual_error, "\n")


# ============================================================
# 9. FINAL SAFE PREDICTION FUNCTION
# ============================================================

predict_with_safety_full <- function(new_expr_scaled,
                                     model_path = "xgb_model_full.rds",
                                     ood_path = "ood_detector_full.rds",
                                     prob_path = "prob_threshold_full.rds",
                                     ent_path = "entropy_threshold_full.rds",
                                     mapping_path = "subtype_name_mapping.rds") {

  model <- readRDS(model_path)
  ood_model <- readRDS(ood_path)
  prob_th <- readRDS(prob_path)
  ent_th  <- readRDS(ent_path)
  name_mapping <- readRDS(mapping_path)

  # OOD detection
  ood_pred <- predict(ood_model, new_expr_scaled)
  in_distribution <- as.logical(ood_pred[1])

  # Prediction (using internally valid class names)
  pred_prob <- predict(model, newdata = new_expr_scaled, type = "prob")
  pred_class_valid <- predict(model, newdata = new_expr_scaled)

  # Restore original names for the output
  pred_class <- name_mapping[as.character(pred_class_valid)]
  colnames(pred_prob) <- name_mapping[colnames(pred_prob)]

  max_prob <- max(pred_prob[1, ])
  entropy  <- -sum(pred_prob[1, ] * log(pred_prob[1, ] + 1e-12))

  # FLAGS (using the original class for messages)
  if (!in_distribution) {
    flag <- "FLAG_OOD"
    final_label <- "EXPERT REVIEW REQUIRED"

  } else if (max_prob < prob_th | entropy > ent_th) {
    flag <- "FLAG_UNCERTAIN"
    final_label <- paste0("POSSIBLE_", pred_class)

  } else {
    flag <- "OK"
    final_label <- as.character(pred_class)
  }

  return(list(
    prediction = final_label,
    raw_class = as.character(pred_class),
    max_probability = round(max_prob, 4),
    entropy = round(entropy, 4),
    flag = flag,
    probabilities = round(pred_prob[1, ], 4)
  ))
}

###############################################################################
####################PREDICTION ON NEW SAMPLE##################################
################################################################################

# ============================================================
# PREDICTION CORRECTED WITH COMPLETE FLAG SYSTEM
# ============================================================

# 1. LOAD ALL MODELS AND PARAMETERS
model <- readRDS("xgb_model_full.rds")
ood_model <- readRDS("ood_detector_full.rds")

scale_means <- readRDS("scale_means.rds")
scale_sds   <- readRDS("scale_sds.rds")

prob_threshold <- readRDS("prob_threshold_full.rds")
entropy_threshold <- readRDS("entropy_threshold_full.rds")

name_mapping <- readRDS("subtype_name_mapping.rds")

# 2. LOAD THE 300 GENES
my_genes <- read.table("Polonen_299_genes_new.txt",
                       header = FALSE,
                       stringsAsFactors = FALSE)[[1]]

# 3. PREPARE NAMES FOR ALIGNMENT
genes_model_Vnames <- setdiff(colnames(model$trainingData), ".outcome")

# Verifications
cat("=== VERIFICATIONS ===\n")
cat("Number of genes in model:", length(genes_model_Vnames), "\n")
cat("Number of genes in Polonen:", length(my_genes), "\n")
cat("Number of genes in scale_means:", length(scale_means), "\n")

# ============================================================
# 4. PROCESS THE NEW SAMPLE (ONCE)
# ============================================================
cat("\n=== PROCESSING SAMPLE ===\n")

new.raw_temp <- read_excel("counts_new.xlsx")
names(new.raw_temp)[1] <- "GeneID"

# Group by GeneID and sum counts (remove duplicates)
new.raw_cleaned <- new.raw_temp %>%
  mutate(across(-GeneID, as.numeric)) %>%
  group_by(GeneID) %>%
  summarise(across(where(is.numeric), sum), .groups = 'drop') %>%
  # Eliminate problematic values
  filter(!is.na(GeneID)) %>%
  filter(GeneID != "") %>%
  filter(!grepl("^\\s*$", GeneID)) %>%
  mutate(GeneID = as.character(GeneID))

# Verifications
cat("Genes after grouping and summing:", nrow(new.raw_cleaned), "\n")
cat("Are there duplicates?", sum(duplicated(new.raw_cleaned$GeneID)) > 0, "\n")

# Convert to matrix with manual row names
gene_names <- new.raw_cleaned$GeneID
new.raw_matrix <- new.raw_cleaned %>%
  select(-GeneID) %>%
  as.matrix()
rownames(new.raw_matrix) <- gene_names

# If there are multiple samples, take the first one
if (ncol(new.raw_matrix) > 1) {
  cat("Multiple samples detected. Using the first column.\n")
  new.raw_matrix <- new.raw_matrix[, 1, drop = FALSE]
}

# Ensure correct orientation (genes in rows)
if (nrow(new.raw_matrix) < ncol(new.raw_matrix)) {
  cat("Transposing matrix to place genes in rows...\n")
  new.raw_matrix <- t(new.raw_matrix)
}

# Clean Ensembl versions
rownames(new.raw_matrix) <- gsub("\\..*", "", rownames(new.raw_matrix))

# ============================================================
# 5. ALIGN WITH THE 300 GENES
# ============================================================
cat("\n=== ALIGNING WITH POLONEN GENES ===\n")

genes_present <- intersect(my_genes, rownames(new.raw_matrix))
genes_missing <- setdiff(my_genes, rownames(new.raw_matrix))

cat("Genes present:", length(genes_present), "\n")
cat("Genes missing:", length(genes_missing), "\n")

if (length(genes_present) == 0) {
  stop("ERROR: NO Polonen genes are present in the sample!")
}

# Extract present genes
expr_sub <- new.raw_matrix[genes_present, , drop = FALSE]

# Add missing genes with value 0
if (length(genes_missing) > 0) {
  fill_mat <- matrix(0,
                     nrow = length(genes_missing),
                     ncol = ncol(expr_sub))
  rownames(fill_mat) <- genes_missing
  colnames(fill_mat) <- colnames(expr_sub)
  expr_sub <- rbind(expr_sub, fill_mat)
}

# Reorder according to my_genes
expr_aligned <- expr_sub[my_genes, , drop = FALSE]

# Check for NAs
if(any(is.na(expr_aligned))) {
  cat("WARNING: NAs found. Replacing with 0.\n")
  expr_aligned[is.na(expr_aligned)] <- 0
}

# Rename rows to V1...V300
rownames(expr_aligned) <- paste0("V", seq_len(nrow(expr_aligned)))

# Transpose to place sample in rows
expression_matrix <- as.data.frame(t(expr_aligned))

# ============================================================
# 6. SCALE DATA
# ============================================================
cat("\n=== SCALING DATA ===\n")

new_expr_scaled <- sweep(expression_matrix, 2, scale_means, "-")
new_expr_scaled <- sweep(new_expr_scaled, 2, scale_sds, "/")
colnames(new_expr_scaled) <- genes_model_Vnames

# Check for NAs after scaling
if(any(is.na(new_expr_scaled))) {
  cat("WARNING: NAs found after scaling. Replacing with 0.\n")
  new_expr_scaled[is.na(new_expr_scaled)] <- 0
}

cat("Final dimensions of scaled data:", dim(new_expr_scaled), "\n")

# ============================================================
# 7. OOD DETECTION
# ============================================================
cat("\n=== OOD DETECTION ===\n")
ood_pred <- predict(ood_model, new_expr_scaled)
in_distribution <- as.logical(ood_pred)
ood_distance <- ifelse(!is.null(attr(ood_pred, "decision.values")),
                       abs(attr(ood_pred, "decision.values")[1]), NA)

cat(ifelse(in_distribution, "[OK]", "[!]"),
    ifelse(in_distribution, "Sample within distribution", "SAMPLE OUT OF DISTRIBUTION (OOD)"), "\n")
if (!is.na(ood_distance)) cat("OOD Distance:", round(ood_distance, 4), "\n")

# ============================================================
# 8. XGBOOST PREDICTION
# ============================================================
cat("\n=== XGBOOST PREDICTION ===\n")
pred_prob <- predict(model, newdata = new_expr_scaled, type = "prob")
pred_class_valid <- predict(model, newdata = new_expr_scaled)

# Apply mapping
pred_class_original <- name_mapping[as.character(pred_class_valid)]
colnames(pred_prob) <- name_mapping[colnames(pred_prob)]

# Metrics
max_prob <- max(pred_prob[1, ])
entropy <- -sum(pred_prob[1, ] * log(pred_prob[1, ] + 1e-12))

cat("Predicted class (original):", pred_class_original, "\n")
cat("Max probability:", round(max_prob, 4), "\n")
cat("Entropy:", round(entropy, 4), "\n")

# ============================================================
# 9. PROBABILITY TABLE
# ============================================================
cat("\n=== PROBABILITIES BY SUBTYPE ===\n")
prob_df <- data.frame(
  Subtype = colnames(pred_prob),
  Probability = round(as.numeric(pred_prob[1,]), 4)
)
prob_df <- prob_df[order(prob_df$Probability, decreasing = TRUE), ]
print(prob_df)

# ============================================================
# 10. HIERARCHICAL FLAG SYSTEM
# ============================================================
cat("\n=== HIERARCHICAL FLAG SYSTEM ===\n")

# Evaluate uncertainty level
uncertainty_level <- case_when(
  max_prob < 0.5 ~ "very high",
  max_prob < 0.7 ~ "high",
  max_prob < prob_threshold ~ "moderate",
  TRUE ~ "low"
)

# HIERARCHICAL FLAGS
if (!in_distribution) {
  # OOD CASE
  if (!is.na(ood_distance)) {
    if (ood_distance > 0.8) {
      flag <- "FLAG_OOD_SEVERE"
      final_label <- "URGENT EXPERT REVIEW REQUIRED"
      recomendacion <- "Highly atypical sample - verify quality"
    } else if (ood_distance > 0.4) {
      flag <- "FLAG_OOD_MODERATE"
      final_label <- paste0("POSSIBLE_", pred_class_original, " (atypical profile)")
      recomendacion <- "Probable batch effect - confirm"
    } else {
      flag <- "FLAG_OOD_MILD"
      final_label <- as.character(pred_class_original)
      recomendacion <- "Slight atypicality - valid with caution"
    }
  } else {
    flag <- "FLAG_OOD"
    final_label <- "EXPERT REVIEW REQUIRED (OOD)"
    recomendacion <- "Sample out of distribution"
  }

} else if (max_prob < prob_threshold | entropy > entropy_threshold) {
  # UNCERTAINTY CASE
  if (max_prob < 0.5) {
    flag <- "FLAG_UNCERTAIN_HIGH"
    final_label <- "UNCLASSIFIABLE - MANDATORY REVIEW"
    recomendacion <- "Model highly insecure (prob < 50%)"
  } else if (entropy > 1.5) {
    flag <- "FLAG_UNCERTAIN_HIGH"
    final_label <- paste0("POSSIBLE_", pred_class_original, " (HIGH UNCERTAINTY)")
    recomendacion <- "High entropy - model confused"
  } else {
    flag <- "FLAG_UNCERTAIN"
    final_label <- paste0("POSSIBLE_", pred_class_original)
    recomendacion <- "Confidence below optimal threshold"
  }

} else {
  # RELIABLE CASE
  flag <- ifelse(max_prob > 0.95, "OK_HIGH_CONFIDENCE", "OK")
  final_label <- as.character(pred_class_original)
  recomendacion <- ifelse(max_prob > 0.95,
                          "High confidence prediction (>95%)",
                          "Reliable prediction")
}

# Display results
cat("\n=== RESULT WITH HIERARCHICAL FLAGS ===\n")
cat("Class:", pred_class_original, "\n")
cat("Max prob:", round(max_prob, 4), "\n")
cat("Entropy:", round(entropy, 4), "\n")
cat("Uncertainty level:", uncertainty_level, "\n")
cat("Flag:", flag, "\n")
cat("Recommendation:", recomendacion, "\n")
cat("Final result:", final_label, "\n")

# Save results
resultado_completo <- list(
  sample_id = "unknown",
  prediction = final_label,
  raw_class = pred_class_original,
  max_probability = round(max_prob, 4),
  entropy = round(entropy, 4),
  flag = flag,
  in_distribution = in_distribution,
  ood_distance = ifelse(!is.na(ood_distance), round(ood_distance, 4), NA),
  uncertainty_level = uncertainty_level,
  recomendacion = recomendacion,
  probabilities = prob_df
)

saveRDS(resultado_completo, paste0("hierarchical_prediction_", Sys.Date(), ".rds"))

# ============================================================
# 11. t-SNE VISUALIZATION
# ============================================================
cat("\n=== GENERATING t-SNE VISUALIZATION ===\n")

# Combine matrices
X_combined <- rbind(expression_matrix_scaled, new_expr_scaled)

# Use original labels
y_combined <- c(as.character(y_original), rep("New", nrow(new_expr_scaled)))

set.seed(42)
tsne_result <- Rtsne(X_combined, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 1000)

tsne_df <- data.frame(
  Dim1 = tsne_result$Y[,1],
  Dim2 = tsne_result$Y[,2],
  Subtype = y_combined
)

# Color palette
colors_subtypes <- c(
  "BCL11B" = "#F8766D",
  "ETP-like" = "#E7851E",
  "HOXA9 TCR" = "#D09400",
  "KMT2A" = "#B2A100",
  "LMO2 gd-like" = "#89AC00",
  "MLLT10" = "#45B500",
  "NKX2-1" = "#00BC51",
  "NKX2-5" = "#00C087",
  "NUP214" = "#00C0B2",
  "NUP98" = "#00C0B3",
  "SPI1" = "#00BCD6",
  "STAG2/LMO2" = "#00B3F2",
  "TAL1 aß-like" = "#28A3FB",
  "TAL1 DP-like" = "#28A3FF",
  "TLX1" = "#9C8DFF",
  "TLX3" = "#D277FF",
  "TME-enriched" = "#F166E8",
  "Others" = "#FF61C7",
  "New" = "red"
)

ggplot(tsne_df, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = Subtype, shape = Subtype), size = 2, alpha = 0.8) +
  scale_color_manual(values = colors_subtypes) +
  scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)) +
  theme_minimal() +
  labs(title = "t-SNE visualization with new samples")

# ============================================================
# 12. PDF PREDICTION TABLE
# ============================================================
cat("\n=== CREATING PDF PROBABILITY TABLE ===\n")

if (nrow(pred_prob) == 1) {
  df <- data.frame(
    Subtype = colnames(pred_prob),
    Probability = round(as.numeric(pred_prob[1, ]), 4)
  )
} else {
  df <- as.data.frame(pred_prob)
  df <- df %>%
    mutate(Subtype = rownames(df)) %>%
    select(Subtype, everything())
  df[, -1] <- round(df[, -1], 4)
}

df <- df %>% arrange(desc(Probability))
print(df)

# Create visual table
tabla_grob <- tableGrob(df, rows = NULL, theme = ttheme_default(
  core = list(fg_params = list(cex = 0.6)),
  colhead = list(fg_params = list(cex = 0.7, fontface = "bold"))
))

ggsave("T_ALL_classifier.pdf", tabla_grob, width = 8, height = 10)
cat("PDF saved: T_ALL_classifier.pdf\n")
