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
# 2. CARGA Y PREPROCESADO
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

# Guardar mapping para restaurar nombres en predicciones
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
# 4. ENTRENAMIENTO XGBOOST
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
# 5. DETECTOR OOD GLOBAL
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
# 6. MÉTRICAS DE INCERTIDUMBRE EN TEST
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
# 7. OPTIMIZACIÓN ROC DE UMBRALES
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

cat("Threshold probabilidad:", prob_threshold, "\n")
cat("Threshold entropía:", entropy_threshold, "\n")

# ============================================================
# 8. EVALUACIÓN COBERTURA vs ERROR RESIDUAL
# ============================================================

prob_threshold <- as.numeric(prob_threshold)
entropy_threshold <- as.numeric(entropy_threshold)

flag_vector <- (uncertainty_df$max_prob < prob_threshold) |
               (uncertainty_df$entropy > entropy_threshold)

coverage <- mean(!flag_vector)
residual_error <- mean(!uncertainty_df$correct[!flag_vector])

cat("Cobertura:", coverage, "\n")
cat("Error residual tras filtro:", residual_error, "\n")


# ============================================================
# 9. FUNCIÓN FINAL DE PREDICCIÓN SEGURA
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
  
  # OOD
  ood_pred <- predict(ood_model, new_expr_scaled)
  in_distribution <- as.logical(ood_pred[1])
  
  # Predicción (con nombres válidos internamente)
  pred_prob <- predict(model, newdata = new_expr_scaled, type = "prob")
  pred_class_valid <- predict(model, newdata = new_expr_scaled)
  
  # Restaurar nombres originales para salida
  pred_class <- name_mapping[as.character(pred_class_valid)]
  colnames(pred_prob) <- name_mapping[colnames(pred_prob)]
  
  max_prob <- max(pred_prob[1, ])
  entropy  <- -sum(pred_prob[1, ] * log(pred_prob[1, ] + 1e-12))
  
  # FLAGS (usando la clase original para mensajes)
  if (!in_distribution) {
    flag <- "FLAG_OOD"
    final_label <- "REVISIÓN EXPERTA REQUERIDA"
    
  } else if (max_prob < prob_th | entropy > ent_th) {
    flag <- "FLAG_UNCERTAIN"
    final_label <- paste0("POSIBLE_", pred_class)
    
  } else {
    flag <- "OK"
    final_label <- as.character(pred_class)
  }
  
  return(list(
    prediction = final_label,
    raw_class = as.character(pred_class),
    max_probability = round(max_prob,4),
    entropy = round(entropy,4),
    flag = flag,
    probabilities = round(pred_prob[1,],4)
  ))
}

###############################################################################
####################PREDICTION ON NEW SAMPLE##################################
################################################################################

# ============================================================
# PREDICCIÓN CORREGIDA CON SISTEMA DE FLAGS COMPLETO
# ============================================================

# 1. CARGAR TODOS LOS MODELOS Y PARÁMETROS
model <- readRDS("xgb_model_full.rds")
ood_model <- readRDS("ood_detector_full.rds")

scale_means <- readRDS("scale_means.rds")
scale_sds   <- readRDS("scale_sds.rds")

prob_threshold <- readRDS("prob_threshold_full.rds")
entropy_threshold <- readRDS("entropy_threshold_full.rds")

name_mapping <- readRDS("subtype_name_mapping.rds")

# 2. CARGAR LOS 300 GENES
my_genes <- read.table("Polonen_299_genes_new.txt", 
                       header = FALSE, 
                       stringsAsFactors = FALSE)[[1]]

# 3. PREPARAR NOMBRES PARA ALINEACIÓN
genes_model_Vnames <- setdiff(colnames(model$trainingData), ".outcome")

# Verificaciones
cat("=== VERIFICACIONES ===\n")
cat("Número de genes en modelo:", length(genes_model_Vnames), "\n")
cat("Número de genes en Polonen:", length(my_genes), "\n")
cat("Número de genes en scale_means:", length(scale_means), "\n")

# ============================================================
# 4. PROCESAR LA NUEVA MUESTRA (UNA SOLA VEZ)
# ============================================================
cat("\n=== PROCESANDO MUESTRA ===\n")

new.raw_temp <- read_excel("counts_new.xlsx")
names(new.raw_temp)[1] <- "GeneID"

# Agrupar por GeneID y sumar conteos (eliminar duplicados)
new.raw_cleaned <- new.raw_temp %>%
  mutate(across(-GeneID, as.numeric)) %>%
  group_by(GeneID) %>%
  summarise(across(where(is.numeric), sum), .groups = 'drop') %>%
  # Eliminar valores problemáticos
  filter(!is.na(GeneID)) %>%
  filter(GeneID != "") %>%
  filter(!grepl("^\\s*$", GeneID)) %>%
  mutate(GeneID = as.character(GeneID))

# Verificaciones
cat("Genes después de agrupar y sumar:", nrow(new.raw_cleaned), "\n")
cat("¿Hay duplicados?", sum(duplicated(new.raw_cleaned$GeneID)) > 0, "\n")

# Convertir a matriz con nombres de fila manualmente
gene_names <- new.raw_cleaned$GeneID
new.raw_matrix <- new.raw_cleaned %>%
  select(-GeneID) %>%
  as.matrix()
rownames(new.raw_matrix) <- gene_names

# Si hay múltiples muestras, tomar la primera
if (ncol(new.raw_matrix) > 1) {
  cat("Múltiples muestras detectadas. Usando la primera columna.\n")
  new.raw_matrix <- new.raw_matrix[, 1, drop = FALSE]
}

# Asegurar orientación correcta (genes en filas)
if (nrow(new.raw_matrix) < ncol(new.raw_matrix)) {
  cat("Transponiendo matriz para tener genes en filas...\n")
  new.raw_matrix <- t(new.raw_matrix)
}

# Limpiar versiones de Ensembl
rownames(new.raw_matrix) <- gsub("\\..*", "", rownames(new.raw_matrix))

# ============================================================
# 5. ALINEAR CON LOS 300 GENES
# ============================================================
cat("\n=== ALINEANDO CON GENES DE POLONEN ===\n")

genes_present <- intersect(my_genes, rownames(new.raw_matrix))
genes_missing <- setdiff(my_genes, rownames(new.raw_matrix))

cat("Genes presentes:", length(genes_present), "\n")
cat("Genes faltantes:", length(genes_missing), "\n")

if (length(genes_present) == 0) {
  stop("¡ERROR: NINGÚN gen de Polonen está presente en la muestra!")
}

# Extraer genes presentes
expr_sub <- new.raw_matrix[genes_present, , drop = FALSE]

# Añadir genes faltantes con valor 0
if (length(genes_missing) > 0) {
  fill_mat <- matrix(0, 
                     nrow = length(genes_missing), 
                     ncol = ncol(expr_sub))
  rownames(fill_mat) <- genes_missing
  colnames(fill_mat) <- colnames(expr_sub)
  expr_sub <- rbind(expr_sub, fill_mat)
}

# Reordenar según my_genes
expr_aligned <- expr_sub[my_genes, , drop = FALSE]

# Verificar NAs
if(any(is.na(expr_aligned))) {
  cat("ADVERTENCIA: NAs encontrados. Reemplazando con 0.\n")
  expr_aligned[is.na(expr_aligned)] <- 0
}

# Renombrar filas a V1...V300
rownames(expr_aligned) <- paste0("V", seq_len(nrow(expr_aligned)))

# Transponer para tener muestra en filas
expression_matrix <- as.data.frame(t(expr_aligned))

# ============================================================
# 6. ESCALAR DATOS
# ============================================================
cat("\n=== ESCALANDO DATOS ===\n")

new_expr_scaled <- sweep(expression_matrix, 2, scale_means, "-")
new_expr_scaled <- sweep(new_expr_scaled, 2, scale_sds, "/")
colnames(new_expr_scaled) <- genes_model_Vnames

# Verificar NAs después del escalado
if(any(is.na(new_expr_scaled))) {
  cat("ADVERTENCIA: NAs después del escalado. Reemplazando con 0.\n")
  new_expr_scaled[is.na(new_expr_scaled)] <- 0
}

cat("Dimensiones finales de datos escalados:", dim(new_expr_scaled), "\n")

# ============================================================
# 7. DETECCIÓN OOD
# ============================================================
cat("\n=== DETECCIÓN OOD ===\n")
ood_pred <- predict(ood_model, new_expr_scaled)
in_distribution <- as.logical(ood_pred)
ood_distance <- ifelse(!is.null(attr(ood_pred, "decision.values")), 
                       abs(attr(ood_pred, "decision.values")[1]), NA)

cat(ifelse(in_distribution, "?", "?"), 
    ifelse(in_distribution, "Muestra dentro de distribución", "MUESTRA FUERA DE DISTRIBUCIÓN (OOD)"), "\n")
if (!is.na(ood_distance)) cat("Distancia OOD:", round(ood_distance, 4), "\n")

# ============================================================
# 8. PREDICCIÓN XGBOOST
# ============================================================
cat("\n=== PREDICCIÓN XGBOOST ===\n")
pred_prob <- predict(model, newdata = new_expr_scaled, type = "prob")
pred_class_valid <- predict(model, newdata = new_expr_scaled)

# Aplicar mapping
pred_class_original <- name_mapping[as.character(pred_class_valid)]
colnames(pred_prob) <- name_mapping[colnames(pred_prob)]

# Métricas
max_prob <- max(pred_prob[1, ])
entropy <- -sum(pred_prob[1, ] * log(pred_prob[1, ] + 1e-12))

cat("Clase predicha (original):", pred_class_original, "\n")
cat("Max probabilidad:", round(max_prob, 4), "\n")
cat("Entropía:", round(entropy, 4), "\n")

# ============================================================
# 9. TABLA DE PROBABILIDADES
# ============================================================
cat("\n=== PROBABILIDADES POR SUBTIPO ===\n")
prob_df <- data.frame(
  Subtype = colnames(pred_prob),
  Probability = round(as.numeric(pred_prob[1,]), 4)
)
prob_df <- prob_df[order(prob_df$Probability, decreasing = TRUE), ]
print(prob_df)

# ============================================================
# 10. SISTEMA DE FLAGS JERÁRQUICOS
# ============================================================
cat("\n=== SISTEMA DE FLAGS JERÁRQUICOS ===\n")

# Evaluar nivel de incertidumbre
uncertainty_level <- case_when(
  max_prob < 0.5 ~ "muy alta",
  max_prob < 0.7 ~ "alta",
  max_prob < prob_threshold ~ "moderada",
  TRUE ~ "baja"
)

# FLAGS JERÁRQUICOS
if (!in_distribution) {
  # CASO OOD
  if (!is.na(ood_distance)) {
    if (ood_distance > 0.8) {
      flag <- "FLAG_OOD_SEVERE"
      final_label <- "REVISIÓN EXPERTA URGENTE"
      recomendacion <- "Muestra muy atípica - verificar calidad"
    } else if (ood_distance > 0.4) {
      flag <- "FLAG_OOD_MODERATE"
      final_label <- paste0("POSIBLE_", pred_class_original, " (perfil atípico)")
      recomendacion <- "Probable batch effect - confirmar"
    } else {
      flag <- "FLAG_OOD_MILD"
      final_label <- as.character(pred_class_original)
      recomendacion <- "Leve atipicidad - válido con precaución"
    }
  } else {
    flag <- "FLAG_OOD"
    final_label <- "REVISIÓN EXPERTA (OOD)"
    recomendacion <- "Muestra fuera de distribución"
  }
  
} else if (max_prob < prob_threshold | entropy > entropy_threshold) {
  # CASO INCERTIDUMBRE
  if (max_prob < 0.5) {
    flag <- "FLAG_UNCERTAIN_HIGH"
    final_label <- "NO CLASIFICABLE - REVISIÓN OBLIGATORIA"
    recomendacion <- "Modelo muy inseguro (prob < 50%)"
  } else if (entropy > 1.5) {
    flag <- "FLAG_UNCERTAIN_HIGH"
    final_label <- paste0("POSIBLE_", pred_class_original, " (ALTA INCERTIDUMBRE)")
    recomendacion <- "Alta entropía - modelo confundido"
  } else {
    flag <- "FLAG_UNCERTAIN"
    final_label <- paste0("POSIBLE_", pred_class_original)
    recomendacion <- "Confianza por debajo del umbral óptimo"
  }
  
} else {
  # CASO CONFIABLE
  flag <- ifelse(max_prob > 0.95, "OK_HIGH_CONFIDENCE", "OK")
  final_label <- as.character(pred_class_original)
  recomendacion <- ifelse(max_prob > 0.95, 
                          "Predicción de alta confianza (>95%)", 
                          "Predicción confiable")
}

# Mostrar resultados
cat("\n=== RESULTADO CON FLAGS JERÁRQUICOS ===\n")
cat("Clase:", pred_class_original, "\n")
cat("Max prob:", round(max_prob, 4), "\n")
cat("Entropía:", round(entropy, 4), "\n")
cat("Nivel incertidumbre:", uncertainty_level, "\n")
cat("Flag:", flag, "\n")
cat("Recomendación:", recomendacion, "\n")
cat("Resultado final:", final_label, "\n")

# Guardar resultado
resultado_completo <- list(
  sample_id = "desconocido",
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

saveRDS(resultado_completo, paste0("prediccion_jerarquica_", Sys.Date(), ".rds"))

# ============================================================
# 11. t-SNE VISUALIZATION
# ============================================================
cat("\n=== GENERANDO VISUALIZACIÓN t-SNE ===\n")

# Combinar matrices
X_combined <- rbind(expression_matrix_scaled, new_expr_scaled)

# Usar etiquetas originales
y_combined <- c(as.character(y_original), rep("New", nrow(new_expr_scaled)))

set.seed(42)
tsne_result <- Rtsne(X_combined, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 1000)

tsne_df <- data.frame(
  Dim1 = tsne_result$Y[,1],
  Dim2 = tsne_result$Y[,2],
  Subtype = y_combined
)

# Paleta de colores
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
# 12. TABLA DE PREDICCIÓN PDF
# ============================================================
cat("\n=== CREANDO TABLA DE PROBABILIDADES PDF ===\n")

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

# Crear tabla gráfica
tabla_grob <- tableGrob(df, rows = NULL, theme = ttheme_default(
  core = list(fg_params = list(cex = 0.6)),
  colhead = list(fg_params = list(cex = 0.7, fontface = "bold"))
))

ggsave("T_ALL_classifier.pdf", tabla_grob, width = 8, height = 10)
cat("PDF guardado: T_ALL_classifier.pdf\n")