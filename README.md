# T-REX
T-REX (T-ALL RNA Expression classifier) is a machine learning-based tool that assigns molecular T-ALL subtypes using bulk RNA sequencing data.  

Project Overview
This tool implements a dual-stream pipeline for the molecular classification of T-cell Acute Lymphoblastic Leukemia (T-ALL)for both predictive accuracy and visual interpretability:

* **Inference Stream: Optimized for speed and direct clinical application using gene-specific scaling.** For classification, raw counts are filtered for the 300-gene panel and then standardized using the training cohort's Z-score parameters (means and SDs). This ensures that the **XGBoost** model receives data in the exact distribution it was trained on, maintaining high sensitivity for subtype identification.

* **Visualization Stream: Utilizes Variance Stabilizing Transformation (VST) to provide high-fidelity clustering and sample-to-cohort comparisons.** For **t-SNE** and **UMAP** projections, data undergoes **Variance Stabilizing Transformation (VST)** via `DESeq2`. This step is crucial for visual analysis as it removes the heteroscedasticity typical of RNA-Seq data (the mean-variance trend), allowing for a clear biological clustering of subtypes without the noise of low-count genes.

By decoupling these processes, the tool provides a high-performance classifier while offering a "gold-standard" visual validation of where each new sample sits within the T-ALL landscape.

Technical Workflow
1. Classification & Inference
The core classifier is built on an XGBoost architecture.

Feature Selection: 300-gene signature proposed by Pölönen et al. 2024.

Preprocessing: New samples undergo automated alignment, missing gene imputation (set to 0), and Z-score normalization using the training cohort's stored means and standard deviations.

Output: Generates a probability distribution across 18 subtypes and a final PDF report.

2. High-Fidelity Visualization
To validate the model's predictions, the pipeline includes a visual QC module:

Normalization: Employs DESeq2::vst to stabilize variance across the transcriptome.

Clustering: Executes t-SNE (T-distributed Stochastic Neighbor Embedding) to project the new sample against the reference cohort, identifying its position relative to known molecular clusters.

Performance & Deliverables
Model: XGBoost trained via caret with 5-fold cross-validation.

Metrics: Full Confusion Matrix and F1-score analysis by subtype are available in the evaluation scripts.

Reporting: Automatic generation of a summary table (T_ALL_classifier.pdf) containing prediction probabilities for each class.

Note: Due to the manuscript being under review at Leukemia (Nature Portfolio), raw data is withheld. The repository includes the trained .rds model and the inference scripts required to process local count matrices.