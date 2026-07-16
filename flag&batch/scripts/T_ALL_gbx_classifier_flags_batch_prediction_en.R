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
# BATCH PREDICTION FUNCTION - CORRECTED VERSION
# ============================================================

predict_batch_safe <- function(new_data_path = "/beegfs_agamede/gscratch/geabialn/T_ALL/Canada_samples.xlsx",
                               model_path = "xgb_model_full.rds",
                               ood_path = "ood_detector_full.rds",
                               prob_path = "prob_threshold_full.rds",
                               ent_path = "entropy_threshold_full.rds",
                               mapping_path = "subtype_name_mapping.rds",
                               genes_path = "Polonen_299_genes_new.txt",
                               output_file = "batch_predictions.csv") {

  # Load models and parameters
  cat("Loading models and parameters...\n")

  # Verify that all required files exist
  required_files <- c(model_path, ood_path, prob_path, ent_path, mapping_path, genes_path)
  for (file in required_files) {
    if (!file.exists(file)) {
      stop(paste("ERROR: File not found:", file))
    }
  }

  model <- readRDS(model_path)
  ood_model <- readRDS(ood_path)

  # Verify scaling files exist
  if (!file.exists("scale_means.rds") | !file.exists("scale_sds.rds")) {
    stop("ERROR: scale_means.rds or scale_sds.rds not found")
  }

  scale_means <- readRDS("scale_means.rds")
  scale_sds <- readRDS("scale_sds.rds")
  prob_th <- readRDS(prob_path)
  ent_th <- readRDS(ent_path)
  name_mapping <- readRDS(mapping_path)
  my_genes <- read.table(genes_path, header = FALSE, stringsAsFactors = FALSE)[[1]]

  # Get gene names from model
  genes_model_Vnames <- setdiff(colnames(model$trainingData), ".outcome")

  # Verify that scaling vector lengths match
  if (length(scale_means) != length(genes_model_Vnames)) {
    cat("WARNING: scale_means length does not match the number of genes\n")
  }

  # Load new data (multiple samples)
  cat("Loading data from:", new_data_path, "\n")

  if (!file.exists(new_data_path)) {
    stop(paste("ERROR: File not found:", new_data_path))
  }

  new.raw_temp <- read_excel(new_data_path)
  new.raw <- as.data.frame(new.raw_temp)
  names(new.raw)[1] <- "GeneID"
  rownames(new.raw) <- new.raw$GeneID
  new.raw <- new.raw[, -1, drop = FALSE]

  # Verify orientation
  cat("Original dimensions:", dim(new.raw), "\n")
  if (ncol(new.raw) > nrow(new.raw)) {
    cat("Transposing matrix...\n")
    new.raw <- as.data.frame(t(new.raw))
  }

  # Clean Ensembl versions
  rownames(new.raw) <- gsub("\\..*", "", rownames(new.raw))
  sample_names <- colnames(new.raw)

  # Align with Polonen genes
  cat("Aligning with reference genes...\n")
  genes_present <- intersect(my_genes, rownames(new.raw))
  genes_missing <- setdiff(my_genes, rownames(new.raw))

  cat("Genes present:", length(genes_present), "/", length(my_genes), "\n")
  cat("Genes missing:", length(genes_missing), "\n")

  if (length(genes_present) == 0) {
    stop("ERROR: NO Polonen genes are present in the samples!")
  }

  # Extract present genes
  expr_sub <- new.raw[genes_present, , drop = FALSE]

  # Add missing genes as zeros
  if (length(genes_missing) > 0) {
    fill_mat <- matrix(0, nrow = length(genes_missing), ncol = ncol(expr_sub))
    rownames(fill_mat) <- genes_missing
    colnames(fill_mat) <- colnames(expr_sub)
    expr_sub <- rbind(as.matrix(expr_sub), fill_mat)
  }

  # Reorder according to my_genes
  expr_aligned <- expr_sub[my_genes, , drop = FALSE]
  rownames(expr_aligned) <- paste0("V", seq_len(nrow(expr_aligned)))

  # Transpose (samples as rows) and convert to data frame
  expression_matrix <- as.data.frame(t(expr_aligned))
  rownames(expression_matrix) <- make.names(sample_names, unique = TRUE)

  # Scale all samples
  cat("Scaling data...\n")

  # Ensure columns match
  if (ncol(expression_matrix) != length(scale_means)) {
    cat("WARNING: Adjusting dimensions for scaling\n")
    # Use only the genes that we have
    common_genes <- intersect(1:ncol(expression_matrix), 1:length(scale_means))
    expression_matrix <- expression_matrix[, common_genes, drop = FALSE]
    scale_means <- scale_means[common_genes]
    scale_sds <- scale_sds[common_genes]
  }

  new_expr_scaled <- sweep(as.matrix(expression_matrix), 2, scale_means, "-")
  new_expr_scaled <- sweep(new_expr_scaled, 2, scale_sds, "/")
  colnames(new_expr_scaled) <- genes_model_Vnames[1:ncol(new_expr_scaled)]

  # Ensure it is a data frame
  new_expr_scaled <- as.data.frame(new_expr_scaled)

  cat("Final matrix dimensions:", dim(new_expr_scaled), "\n")
  cat("Number of samples to predict:", nrow(new_expr_scaled), "\n")

  # ============================================================
  # BATCH PREDICTIONS
  # ============================================================

  cat("\n=== STARTING BATCH PREDICTIONS ===\n")

  # 1. OOD detection for all samples
  cat("Running OOD detection...\n")
  ood_pred_all <- tryCatch({
    predict(ood_model, new_expr_scaled)
  }, error = function(e) {
    cat("Error in OOD prediction:", e$message, "\n")
    cat("Using fallback method...\n")
    # Fallback: flag all as true (in-distribution)
    rep(TRUE, nrow(new_expr_scaled))
  })

  in_distribution_all <- as.logical(ood_pred_all)

  # Get OOD distances if available
  ood_distances <- if (!is.null(attr(ood_pred_all, "decision.values"))) {
    as.numeric(abs(attr(ood_pred_all, "decision.values")[, 1]))
  } else {
    rep(NA, nrow(new_expr_scaled))
  }

  # 2. XGBoost Predictions (probabilities and classes)
  cat("Running XGBoost predictions...\n")
  pred_prob_all <- tryCatch({
    predict(model, newdata = new_expr_scaled, type = "prob")
  }, error = function(e) {
    cat("Error in probability prediction:", e$message, "\n")
    # Create dummy probability matrix
    n_classes <- length(unique(model$finalModel$class_names))
    matrix(1/n_classes, nrow = nrow(new_expr_scaled), ncol = n_classes)
  })

  pred_class_valid_all <- predict(model, newdata = new_expr_scaled)

  # Ensure pred_prob_all is a data frame
  pred_prob_all <- as.data.frame(pred_prob_all)

  # 3. Calculate uncertainty metrics
  cat("Calculating uncertainty metrics...\n")
  max_prob_all <- apply(pred_prob_all, 1, max, na.rm = TRUE)
  entropy_all <- apply(pred_prob_all, 1, function(p) {
    p <- as.numeric(p)
    p <- p[!is.na(p)]
    if (length(p) == 0) return(NA)
    -sum(p * log(p + 1e-12), na.rm = TRUE)
  })

  # 4. Map back to original names
  pred_class_original_all <- name_mapping[as.character(pred_class_valid_all)]

  # Rename probability columns
  if (!is.null(colnames(pred_prob_all))) {
    colnames(pred_prob_all) <- name_mapping[colnames(pred_prob_all)]
  }

  # ============================================================
  # HIERARCHICAL FLAG SYSTEM (VECTORIZED)
  # ============================================================

  cat("\n=== APPLYING FLAG SYSTEM ===\n")

  n_samples <- nrow(new_expr_scaled)

  # Initialize results vectors
  flag_vector <- character(n_samples)
  final_label_vector <- character(n_samples)
  recommendation_vector <- character(n_samples)
  uncertainty_level_vector <- character(n_samples)

  # Evaluate uncertainty levels for all samples
  for (i in 1:n_samples) {
    if (is.na(max_prob_all[i])) {
      uncertainty_level_vector[i] <- "unknown"
    } else if (max_prob_all[i] < 0.5) {
      uncertainty_level_vector[i] <- "very high"
    } else if (max_prob_all[i] < 0.7) {
      uncertainty_level_vector[i] <- "high"
    } else if (max_prob_all[i] < prob_th) {
      uncertainty_level_vector[i] <- "moderate"
    } else {
      uncertainty_level_vector[i] <- "low"
    }
  }

  # Apply flagging rules
  for (i in 1:n_samples) {
    in_dist <- ifelse(is.na(in_distribution_all[i]), FALSE, in_distribution_all[i])
    max_prob <- ifelse(is.na(max_prob_all[i]), 0, max_prob_all[i])
    entropy <- ifelse(is.na(entropy_all[i]), Inf, entropy_all[i])
    ood_dist <- ifelse(is.na(ood_distances[i]), 0, ood_distances[i])
    pred_class <- ifelse(is.na(pred_class_original_all[i]), "UNKNOWN", pred_class_original_all[i])

    if (!in_dist) {
      # OOD CASE
      if (!is.na(ood_dist)) {
        if (ood_dist > 0.8) {
          flag_vector[i] <- "FLAG_OOD_SEVERE"
          final_label_vector[i] <- "URGENT EXPERT REVIEW"
          recommendation_vector[i] <- "Highly atypical sample - verify quality"
        } else if (ood_dist > 0.4) {
          flag_vector[i] <- "FLAG_OOD_MODERATE"
          final_label_vector[i] <- paste0("POSSIBLE_", pred_class, " (atypical profile)")
          recommendation_vector[i] <- "Probable batch effect - confirm"
        } else {
          flag_vector[i] <- "FLAG_OOD_MILD"
          final_label_vector[i] <- pred_class
          recommendation_vector[i] <- "Slight atypicality - valid with caution"
        }
      } else {
        flag_vector[i] <- "FLAG_OOD"
        final_label_vector[i] <- "EXPERT REVIEW (OOD)"
        recommendation_vector[i] <- "Sample out of distribution"
      }

    } else if (max_prob < prob_th | entropy > ent_th) {
      # UNCERTAINTY CASE
      if (max_prob < 0.5) {
        flag_vector[i] <- "FLAG_UNCERTAIN_HIGH"
        final_label_vector[i] <- "UNCLASSIFIABLE - MANDATORY REVIEW"
        recommendation_vector[i] <- "Model highly insecure (prob < 50%)"
      } else if (entropy > 1.5) {
        flag_vector[i] <- "FLAG_UNCERTAIN_HIGH"
        final_label_vector[i] <- paste0("POSSIBLE_", pred_class, " (HIGH UNCERTAINTY)")
        recommendation_vector[i] <- "High entropy - model confused"
      } else {
        flag_vector[i] <- "FLAG_UNCERTAIN"
        final_label_vector[i] <- paste0("POSSIBLE_", pred_class)
        recommendation_vector[i] <- "Confidence below optimal threshold"
      }

    } else {
      # RELIABLE CASE
      flag_vector[i] <- ifelse(max_prob > 0.95, "OK_HIGH_CONFIDENCE", "OK")
      final_label_vector[i] <- pred_class
      recommendation_vector[i] <- ifelse(max_prob > 0.95,
                                         "High confidence prediction (>95%)",
                                         "Reliable prediction")
    }
  }

  # ============================================================
  # CREATE RESULTS TABLE
  # ============================================================

  cat("\n=== CREATING RESULTS TABLE ===\n")

  # Summary Data Frame
  resultados_df <- data.frame(
    Sample_ID = sample_names,
    Predicted_Class = pred_class_original_all,
    Final_Label = final_label_vector,
    Max_Probability = round(as.numeric(max_prob_all), 4),
    Entropy = round(as.numeric(entropy_all), 4),
    In_Distribution = in_distribution_all,
    OOD_Distance = round(as.numeric(ood_distances), 4),
    Uncertainty_Level = uncertainty_level_vector,
    Flag = flag_vector,
    Recommendation = recomendacion_vector,
    stringsAsFactors = FALSE
  )

  # Append top 3 probabilities for each sample
  cat("Appending top probabilities per sample...\n")
  top_probs_list <- character(n_samples)
  for (i in 1:n_samples) {
    if (i <= nrow(pred_prob_all)) {
      sample_probs <- as.numeric(pred_prob_all[i, ])
      if (length(sample_probs) > 0 && !all(is.na(sample_probs))) {
        names(sample_probs) <- colnames(pred_prob_all)
        top3_indices <- order(sample_probs, decreasing = TRUE)[1:min(3, length(sample_probs))]
        top3_values <- sample_probs[top3_indices]
        top3_names <- names(top3_values)
        top3_str <- paste0(paste0(top3_names, " (", round(top3_values, 3), ")"), collapse = "; ")
        top_probs_list[i] <- top3_str
      } else {
        top_probs_list[i] <- "Not available"
      }
    } else {
      top_probs_list[i] <- "Not available"
    }
  }
  resultados_df$Top3_Probabilities <- top_probs_list

  # Complete probability matrix
  if (nrow(pred_prob_all) > 0) {
    prob_matrix <- as.data.frame(pred_prob_all)
    prob_matrix$Sample_ID <- sample_names[1:nrow(prob_matrix)]
    prob_matrix <- prob_matrix[, c("Sample_ID", setdiff(colnames(prob_matrix), "Sample_ID"))]
  } else {
    prob_matrix <- data.frame(Sample_ID = character(0))
  }

  # ============================================================
  # SAVE RESULTS
  # ============================================================

  # Save summary CSV
  write.csv(resultados_df, output_file, row.names = FALSE)
  cat("Results saved to:", output_file, "\n")

  # Save full probability matrix
  if (nrow(prob_matrix) > 0) {
    prob_file <- gsub(".csv", "_probabilities.csv", output_file)
    write.csv(prob_matrix, prob_file, row.names = FALSE)
    cat("Probability matrix saved to:", prob_file, "\n")
  }

  # Save complete R object for downstream analysis
  rds_file <- gsub(".csv", ".rds", output_file)
  complete_results <- list(
    summary = resultados_df,
    probabilities = prob_matrix,
    raw_predictions = if(exists("pred_class_valid_all")) pred_class_valid_all else NULL,
    raw_probabilities = pred_prob_all,
    metadata = list(
      n_samples = n_samples,
      date = Sys.Date(),
      model_used = model_path,
      thresholds = c(prob_threshold = prob_th, entropy_threshold = ent_th)
    )
  )
  saveRDS(complete_results, rds_file)
  cat("Complete R object saved to:", rds_file, "\n")

  # ============================================================
  # VISUALIZATIONS
  # ============================================================

  cat("\n=== GENERATING VISUALIZATIONS ===\n")

  tryCatch({
    # Confidence histogram
    p1 <- ggplot(resultados_df, aes(x = Max_Probability, fill = In_Distribution)) +
      geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
      geom_vline(xintercept = prob_th, linetype = "dashed", color = "red") +
      theme_minimal() +
      labs(title = "Distribution of Maximum Probabilities",
           x = "Maximum Probability", y = "Frequency")

    # Flag distribution
    flag_counts <- as.data.frame(table(Flag = resultados_df$Flag))
    if (nrow(flag_counts) > 0) {
      p2 <- ggplot(flag_counts, aes(x = reorder(Flag, -Freq), y = Freq, fill = Flag)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Flag Distribution",
             x = "Flag", y = "Number of Samples")

      # Save visualizations
      pdf_file <- gsub(".csv", "_visualizations.pdf", output_file)
      pdf(pdf_file, width = 10, height = 8)
      grid.arrange(p1, p2, ncol = 1)
      dev.off()
      cat("Visualizations saved to:", pdf_file, "\n")
    }
  }, error = function(e) {
    cat("Could not generate visualizations:", e$message, "\n")
  })

  # Display summary
  cat("\n=== PREDICTION SUMMARY ===\n")
  cat("Total processed samples:", n_samples, "\n")
  cat("\nFlag distribution:\n")
  print(table(flag_vector))
  cat("\nFirst 5 samples:\n")
  print(head(resultados_df[, c("Sample_ID", "Final_Label", "Flag", "Max_Probability")], 5))

  return(resultados_df)
}

# ============================================================
# EXECUTE FUNCTION
# ============================================================

# Set working directory
setwd("/beegfs_agamede/gscratch/geabialn/T_ALL/")

# Run batch predictions
resultados <- predict_batch_safe()
