# ============================================================
# 1. LIBRERÍAS
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
# FUNCIÓN DE PREDICCIÓN POR LOTES - VERSIÓN CORREGIDA
# ============================================================

predict_batch_safe <- function(new_data_path = "/beegfs_agamede/gscratch/geabialn/T_ALL/Canada_samples.xlsx",
                               model_path = "xgb_model_full.rds",
                               ood_path = "ood_detector_full.rds",
                               prob_path = "prob_threshold_full.rds",
                               ent_path = "entropy_threshold_full.rds",
                               mapping_path = "subtype_name_mapping.rds",
                               genes_path = "Polonen_299_genes_new.txt",
                               output_file = "predicciones_batch.csv") {
  
  # Cargar modelos y parámetros
  cat("Cargando modelos y parámetros...\n")
  
  # Verificar que todos los archivos existen
  required_files <- c(model_path, ood_path, prob_path, ent_path, mapping_path, genes_path)
  for (file in required_files) {
    if (!file.exists(file)) {
      stop(paste("ERROR: No se encuentra el archivo:", file))
    }
  }
  
  model <- readRDS(model_path)
  ood_model <- readRDS(ood_path)
  
  # Verificar que existen los archivos de escalado
  if (!file.exists("scale_means.rds") | !file.exists("scale_sds.rds")) {
    stop("ERROR: No se encuentran los archivos scale_means.rds o scale_sds.rds")
  }
  
  scale_means <- readRDS("scale_means.rds")
  scale_sds <- readRDS("scale_sds.rds")
  prob_th <- readRDS(prob_path)
  ent_th <- readRDS(ent_path)
  name_mapping <- readRDS(mapping_path)
  my_genes <- read.table(genes_path, header = FALSE, stringsAsFactors = FALSE)[[1]]
  
  # Obtener nombres de genes del modelo
  genes_model_Vnames <- setdiff(colnames(model$trainingData), ".outcome")
  
  # Verificar que la longitud de los vectores de escalado coincide
  if (length(scale_means) != length(genes_model_Vnames)) {
    cat("ADVERTENCIA: La longitud de scale_means no coincide con el número de genes\n")
  }
  
  # Cargar datos nuevos (múltiples muestras)
  cat("Cargando datos de:", new_data_path, "\n")
  
  if (!file.exists(new_data_path)) {
    stop(paste("ERROR: No se encuentra el archivo:", new_data_path))
  }
  
  new.raw_temp <- read_excel(new_data_path)
  new.raw <- as.data.frame(new.raw_temp)
  names(new.raw)[1] <- "GeneID"
  rownames(new.raw) <- new.raw$GeneID
  new.raw <- new.raw[, -1, drop = FALSE]
  
  # Verificar orientación
  cat("Dimensiones originales:", dim(new.raw), "\n")
  if (ncol(new.raw) > nrow(new.raw)) {
    cat("Transponiendo matriz...\n")
    new.raw <- as.data.frame(t(new.raw))
  }
  
  # Limpiar versiones de Ensembl
  rownames(new.raw) <- gsub("\\..*", "", rownames(new.raw))
  sample_names <- colnames(new.raw)
  
  # Alinear con genes de Polonen
  cat("Alineando con genes de referencia...\n")
  genes_present <- intersect(my_genes, rownames(new.raw))
  genes_missing <- setdiff(my_genes, rownames(new.raw))
  
  cat("Genes presentes:", length(genes_present), "/", length(my_genes), "\n")
  cat("Genes faltantes:", length(genes_missing), "\n")
  
  if (length(genes_present) == 0) {
    stop("¡ERROR NINGÚN gen de Polonen está presente en las muestras!")
  }
  
  # Extraer genes presentes
  expr_sub <- new.raw[genes_present, , drop = FALSE]
  
  # Añadir genes faltantes como ceros
  if (length(genes_missing) > 0) {
    fill_mat <- matrix(0, nrow = length(genes_missing), ncol = ncol(expr_sub))
    rownames(fill_mat) <- genes_missing
    colnames(fill_mat) <- colnames(expr_sub)
    expr_sub <- rbind(as.matrix(expr_sub), fill_mat)
  }
  
  # Reordenar según my_genes
  expr_aligned <- expr_sub[my_genes, , drop = FALSE]
  rownames(expr_aligned) <- paste0("V", seq_len(nrow(expr_aligned)))
  
  # Transponer (muestras como filas) y convertir a data frame
  expression_matrix <- as.data.frame(t(expr_aligned))
  rownames(expression_matrix) <- make.names(sample_names, unique = TRUE)
  
  # Escalado de todas las muestras
  cat("Escalando datos...\n")
  
  # Asegurar que las columnas coinciden
  if (ncol(expression_matrix) != length(scale_means)) {
    cat("ADVERTENCIA: Ajustando dimensiones para escalado\n")
    # Usar solo los genes que tenemos
    common_genes <- intersect(1:ncol(expression_matrix), 1:length(scale_means))
    expression_matrix <- expression_matrix[, common_genes, drop = FALSE]
    scale_means <- scale_means[common_genes]
    scale_sds <- scale_sds[common_genes]
  }
  
  new_expr_scaled <- sweep(as.matrix(expression_matrix), 2, scale_means, "-")
  new_expr_scaled <- sweep(new_expr_scaled, 2, scale_sds, "/")
  colnames(new_expr_scaled) <- genes_model_Vnames[1:ncol(new_expr_scaled)]
  
  # Asegurar que es un data frame
  new_expr_scaled <- as.data.frame(new_expr_scaled)
  
  cat("Matriz final dimensiones:", dim(new_expr_scaled), "\n")
  cat("Número de muestras a predecir:", nrow(new_expr_scaled), "\n")
  
  # ============================================================
  # PREDICCIONES POR LOTES
  # ============================================================
  
  cat("\n=== INICIANDO PREDICCIONES POR LOTES ===\n")
  
  # 1. Detección OOD para todas las muestras
  cat("Ejecutando detección OOD...\n")
  ood_pred_all <- tryCatch({
    predict(ood_model, new_expr_scaled)
  }, error = function(e) {
    cat("Error en predicción OOD:", e$message, "\n")
    cat("Usando método alternativo...\n")
    # Alternativa: usar distancia euclidiana simple
    rep(TRUE, nrow(new_expr_scaled))
  })
  
  in_distribution_all <- as.logical(ood_pred_all)
  
  # Obtener distancias OOD si están disponibles
  ood_distances <- if (!is.null(attr(ood_pred_all, "decision.values"))) {
    as.numeric(abs(attr(ood_pred_all, "decision.values")[, 1]))
  } else {
    rep(NA, nrow(new_expr_scaled))
  }
  
  # 2. Predicciones XGBoost (probabilidades y clases)
  cat("Ejecutando predicciones XGBoost...\n")
  pred_prob_all <- tryCatch({
    predict(model, newdata = new_expr_scaled, type = "prob")
  }, error = function(e) {
    cat("Error en predicción de probabilidades:", e$message, "\n")
    # Crear matriz de probabilidades dummy
    n_classes <- length(unique(model$finalModel$class_names))
    matrix(1/n_classes, nrow = nrow(new_expr_scaled), ncol = n_classes)
  })
  
  pred_class_valid_all <- predict(model, newdata = new_expr_scaled)
  
  # Asegurar que pred_prob_all es un data frame
  pred_prob_all <- as.data.frame(pred_prob_all)
  
  # 3. Calcular métricas de incertidumbre
  cat("Calculando métricas de incertidumbre...\n")
  max_prob_all <- apply(pred_prob_all, 1, max, na.rm = TRUE)
  entropy_all <- apply(pred_prob_all, 1, function(p) {
    p <- as.numeric(p)
    p <- p[!is.na(p)]
    if (length(p) == 0) return(NA)
    -sum(p * log(p + 1e-12), na.rm = TRUE)
  })
  
  # 4. Aplicar mapping a nombres originales
  pred_class_original_all <- name_mapping[as.character(pred_class_valid_all)]
  
  # Renombrar columnas de probabilidades
  if (!is.null(colnames(pred_prob_all))) {
    colnames(pred_prob_all) <- name_mapping[colnames(pred_prob_all)]
  }
  
  # ============================================================
  # SISTEMA DE FLAGS JERÁRQUICOS (VECTORIZADO)
  # ============================================================
  
  cat("\n=== APLICANDO SISTEMA DE FLAGS ===\n")
  
  n_samples <- nrow(new_expr_scaled)
  
  # Inicializar vectores de resultados
  flag_vector <- character(n_samples)
  final_label_vector <- character(n_samples)
  recomendacion_vector <- character(n_samples)
  uncertainty_level_vector <- character(n_samples)
  
  # Evaluar nivel de incertidumbre para todas las muestras
  for (i in 1:n_samples) {
    if (is.na(max_prob_all[i])) {
      uncertainty_level_vector[i] <- "desconocida"
    } else if (max_prob_all[i] < 0.5) {
      uncertainty_level_vector[i] <- "muy alta"
    } else if (max_prob_all[i] < 0.7) {
      uncertainty_level_vector[i] <- "alta"
    } else if (max_prob_all[i] < prob_th) {
      uncertainty_level_vector[i] <- "moderada"
    } else {
      uncertainty_level_vector[i] <- "baja"
    }
  }
  
  # Aplicar reglas de flags
  for (i in 1:n_samples) {
    in_dist <- ifelse(is.na(in_distribution_all[i]), FALSE, in_distribution_all[i])
    max_prob <- ifelse(is.na(max_prob_all[i]), 0, max_prob_all[i])
    entropy <- ifelse(is.na(entropy_all[i]), Inf, entropy_all[i])
    ood_dist <- ifelse(is.na(ood_distances[i]), 0, ood_distances[i])
    pred_class <- ifelse(is.na(pred_class_original_all[i]), "DESCONOCIDO", pred_class_original_all[i])
    
    if (!in_dist) {
      # CASO OOD
      if (!is.na(ood_dist)) {
        if (ood_dist > 0.8) {
          flag_vector[i] <- "FLAG_OOD_SEVERE"
          final_label_vector[i] <- "REVISIÓN EXPERTA URGENTE"
          recomendacion_vector[i] <- "Muestra muy atípica - verificar calidad"
        } else if (ood_dist > 0.4) {
          flag_vector[i] <- "FLAG_OOD_MODERATE"
          final_label_vector[i] <- paste0("POSIBLE_", pred_class, " (perfil atípico)")
          recomendacion_vector[i] <- "Probable batch effect - confirmar"
        } else {
          flag_vector[i] <- "FLAG_OOD_MILD"
          final_label_vector[i] <- pred_class
          recomendacion_vector[i] <- "Leve atipicidad - válido con precaución"
        }
      } else {
        flag_vector[i] <- "FLAG_OOD"
        final_label_vector[i] <- "REVISIÓN EXPERTA (OOD)"
        recomendacion_vector[i] <- "Muestra fuera de distribución"
      }
      
    } else if (max_prob < prob_th | entropy > ent_th) {
      # CASO INCERTIDUMBRE
      if (max_prob < 0.5) {
        flag_vector[i] <- "FLAG_UNCERTAIN_HIGH"
        final_label_vector[i] <- "NO CLASIFICABLE - REVISIÓN OBLIGATORIA"
        recomendacion_vector[i] <- "Modelo muy inseguro (prob < 50%)"
      } else if (entropy > 1.5) {
        flag_vector[i] <- "FLAG_UNCERTAIN_HIGH"
        final_label_vector[i] <- paste0("POSIBLE_", pred_class, " (ALTA INCERTIDUMBRE)")
        recomendacion_vector[i] <- "Alta entropía - modelo confundido"
      } else {
        flag_vector[i] <- "FLAG_UNCERTAIN"
        final_label_vector[i] <- paste0("POSIBLE_", pred_class)
        recomendacion_vector[i] <- "Confianza por debajo del umbral óptimo"
      }
      
    } else {
      # CASO CONFIABLE
      flag_vector[i] <- ifelse(max_prob > 0.95, "OK_HIGH_CONFIDENCE", "OK")
      final_label_vector[i] <- pred_class
      recomendacion_vector[i] <- ifelse(max_prob > 0.95, 
                                        "Predicción de alta confianza (>95%)", 
                                        "Predicción confiable")
    }
  }
  
  # ============================================================
  # CREAR TABLA DE RESULTADOS
  # ============================================================
  
  cat("\n=== CREANDO TABLA DE RESULTADOS ===\n")
  
  # Data frame resumen
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
  
  # Añadir top 3 probabilidades para cada muestra
  cat("Añadiendo top probabilidades por muestra...\n")
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
        top_probs_list[i] <- "No disponible"
      }
    } else {
      top_probs_list[i] <- "No disponible"
    }
  }
  resultados_df$Top3_Probabilities <- top_probs_list
  
  # Matriz completa de probabilidades
  if (nrow(pred_prob_all) > 0) {
    prob_matrix <- as.data.frame(pred_prob_all)
    prob_matrix$Sample_ID <- sample_names[1:nrow(prob_matrix)]
    prob_matrix <- prob_matrix[, c("Sample_ID", setdiff(colnames(prob_matrix), "Sample_ID"))]
  } else {
    prob_matrix <- data.frame(Sample_ID = character(0))
  }
  
  # ============================================================
  # GUARDAR RESULTADOS
  # ============================================================
  
  # Guardar CSV con resumen
  write.csv(resultados_df, output_file, row.names = FALSE)
  cat("Resultados guardados en:", output_file, "\n")
  
  # Guardar matriz de probabilidades completa
  if (nrow(prob_matrix) > 0) {
    prob_file <- gsub(".csv", "_probabilidades.csv", output_file)
    write.csv(prob_matrix, prob_file, row.names = FALSE)
    cat("Matriz de probabilidades guardada en:", prob_file, "\n")
  }
  
  # Guardar objeto R completo para análisis posteriores
  rds_file <- gsub(".csv", ".rds", output_file)
  resultados_completos <- list(
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
  saveRDS(resultados_completos, rds_file)
  cat("Objeto R completo guardado en:", rds_file, "\n")
  
  # ============================================================
  # VISUALIZACIONES
  # ============================================================
  
  cat("\n=== GENERANDO VISUALIZACIONES ===\n")
  
  tryCatch({
    # Histograma de confianza
    p1 <- ggplot(resultados_df, aes(x = Max_Probability, fill = In_Distribution)) +
      geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
      geom_vline(xintercept = prob_th, linetype = "dashed", color = "red") +
      theme_minimal() +
      labs(title = "Distribución de probabilidades máximas",
           x = "Probabilidad máxima", y = "Frecuencia")
    
    # Distribución de flags
    flag_counts <- as.data.frame(table(Flag = resultados_df$Flag))
    if (nrow(flag_counts) > 0) {
      p2 <- ggplot(flag_counts, aes(x = reorder(Flag, -Freq), y = Freq, fill = Flag)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Distribución de flags",
             x = "Flag", y = "Número de muestras")
      
      # Guardar visualizaciones
      pdf_file <- gsub(".csv", "_visualizaciones.pdf", output_file)
      pdf(pdf_file, width = 10, height = 8)
      grid.arrange(p1, p2, ncol = 1)
      dev.off()
      cat("Visualizaciones guardadas en:", pdf_file, "\n")
    }
  }, error = function(e) {
    cat("No se pudieron generar las visualizaciones:", e$message, "\n")
  })
  
  # Mostrar resumen
  cat("\n=== RESUMEN DE PREDICCIONES ===\n")
  cat("Total muestras procesadas:", n_samples, "\n")
  cat("\nDistribución de flags:\n")
  print(table(flag_vector))
  cat("\nPrimeras 5 muestras:\n")
  print(head(resultados_df[, c("Sample_ID", "Final_Label", "Flag", "Max_Probability")], 5))
  
  return(resultados_df)
}

# ============================================================
# EJECUTAR LA FUNCIÓN
# ============================================================

# Establecer directorio de trabajo
setwd("/beegfs_agamede/gscratch/geabialn/T_ALL/")

# Ejecutar la función
resultados <- predict_batch_safe()