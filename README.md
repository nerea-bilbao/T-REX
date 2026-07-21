# T-REX (T-ALL RNA Expression Classifier)

**T-REX** is a machine learning-based tool designed to assign molecular **T-cell Acute Lymphoblastic Leukemia (T-ALL)** subtypes using bulk RNA sequencing data.

---

# Project Overview

T-REX implements a dual-stream pipeline for the molecular classification of T-ALL, balancing predictive accuracy with high-fidelity visual interpretability.

## Inference Stream

Optimized for speed and direct clinical application using gene-specific scaling.

For classification:

- Raw counts are filtered for the **300-gene panel**.
- Gene expression is standardized using the training cohort's **Z-score parameters** (means and standard deviations).
- The processed data are then supplied to the **XGBoost** classifier using the same feature distribution employed during model training.

This approach maintains high sensitivity for molecular subtype identification.

## Visualization Stream

The visualization workflow uses **t-distributed Stochastic Neighbor Embedding (t-SNE)** on normalized expression data to project new samples alongside the reference cohort.

This provides:

- High-fidelity visual clustering.
- Intuitive comparison between new samples and reference cases.
- Visual validation of subtype assignment.

By decoupling prediction and visualization, T-REX provides both a high-performance classifier and a "gold-standard" visual representation of each sample within the T-ALL molecular landscape.

---

# Running Modes

## Single-Sample Mode (`/single`)

### Purpose

Fast classification of individual samples.

### Best for

- Real-time analysis.
- Single-patient clinical cases.

### Workflow

Processes count matrices one sample at a time to produce:

- Molecular subtype prediction.
- Quality-control visualization.

---

## Batch & Uncertainty Mode (`/flag&batch`)

### Purpose

Cohort-level classification with prediction reliability estimation.

### Best for

- Research cohorts.
- Large datasets.
- Analyses where confidence estimation is important.

### Workflow

Processes multiple samples simultaneously and computes uncertainty metrics to:

- Flag ambiguous predictions.
- Detect potential outliers.
- Identify samples requiring additional validation.

---

# Input Data Requirements

Before running the pipeline, ensure that the appropriate `input/` directory (inside either `single/` or `flag&batch/`) contains the following files.

| File | Description |
|------|-------------|
| `counts_new.xlsx` | User count matrix containing the samples to classify. |
| `Polonen_299_genes_new.txt` | Target gene signature list. |
| `T_ALL_labels_full_cohort_total_ordered.xlsx` | Reference cohort subtype annotations and metadata. |
| `data_counts_new` | Reference cohort count matrix (download required). |

---

# Required External Download

The reference expression dataset (`data_counts_new`) is hosted on **Synapse** due to file size limitations.

Please download the bulk RNA-seq count matrix from the official repository:

> **Pölönen et al. Synapse Repository**  
> Synapse ID: **syn54032669**

---

# Technical Workflow

## Classification & Inference

### Machine Learning Model

- **Algorithm:** XGBoost
- **Training framework:** `caret`
- **Validation:** 5-fold cross-validation

### Feature Selection

- 300-gene molecular signature proposed by **Pölönen et al. (2024)**.

### Preprocessing

The pipeline performs:

1. Automatic gene alignment.
2. Missing gene imputation (missing genes are assigned a value of 0).
3. Z-score normalization using stored training parameters.

### Output

The classifier returns a probability distribution across **18 molecular T-ALL subtypes**.

---

## High-Fidelity Visualization

### Clustering

Samples are projected against the reference cohort using **t-SNE**, allowing visual assessment of subtype similarity and sample positioning.

---

# Deliverables

## Summary Report

An automatic PDF report named:

```text
T_ALL_classifier.pdf
```

is generated and includes:

- Predicted subtype.
- Probability for each molecular subtype.
- Visualization of the sample within the reference cohort.

## Model Files

Pre-trained `.rds` models and inference scripts are included within the corresponding execution directories:

- `single/`
- `flag&batch/`

---

# Repository Structure

```text
T-REX/
├── single/          # Single-sample processing mode
└── flag&batch/      # Batch processing mode with uncertainty metrics
```

Both directories contain the necessary scripts and resources to process local RNA-seq count matrices.

---

# References

- Pölönen et al. (2024). Molecular classification of T-cell Acute Lymphoblastic Leukemia.
- Synapse Repository: **syn54032669**

