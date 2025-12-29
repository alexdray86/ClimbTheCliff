# ClimbTheCliff

An R package for deconvolution of bulk RNA-seq data using single-cell RNA-seq references, with specialized methods for cell-type abundance estimation, high-resolution cell-type-specific gene expression prediction, and drug sensitivity analysis.

## Installation

You can install the development version of ClimbTheCliff from GitHub using either `devtools` or `remotes`:

### Using devtools

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install ClimbTheCliff
devtools::install_github("drayac/ClimbTheCliff", force=T)
```

### Using remotes (alternative)

```r
# Install remotes if you haven't already
install.packages("remotes")

# Install ClimbTheCliff
remotes::install_github("drayac/ClimbTheCliff", force=T)
```

## Load the package

```r
library(ClimbTheCliff)
```

---

## Methods Overview

The ClimbTheCliff package provides two main deconvolution methods:

1. **CLIMB** - Cell-type abundance and expression deconvolution
2. **CLIFF** - Drug sensitivity analysis using deconvoluted cell-type expression

---

## CLIMB: Cell-type Abundance and Expression Deconvolution

### Description

CLIMB (Ce**L**l-type abundance and expression deconvolution using **I**nverse **M**odeling of **B**ulk RNA-seq) is a method that deconvolutes bulk RNA-seq data using a single-cell RNA-seq reference dataset. It can predict:

- **Cell-type proportions** in bulk samples
- **High-resolution cell-type-specific gene expression** for each sample
- **Differential expression** between conditions or at the sample level

### Key Features

- Uses an empirical Bayes procedure with subsampling for robust estimation
- Handles cancer-specific differential expression
- Supports multiple operation modes for different analysis needs
- Provides both overall and sample-specific cell-type expression profiles

### Usage

```r
climb(sc, bulk, mode = "abundance", up.lim = Inf, lambda = 0, 
      verbose = TRUE, cancer_pattern = "*", conditions = NA, 
      final_res = list(), min_common_genes = 100, 
      ratio_cell_increase = 0.02, n.iter.subsampling = 5, 
      min.n.cells = 50, n.top_mean.genes = 500)
```

### Parameters

- **`sc`**: ExpressionSet object containing the single-cell RNA-seq reference dataset
- **`bulk`**: ExpressionSet object containing the bulk RNA-seq mixtures to be deconvoluted
- **`mode`**: Operating mode (see below for details)
  - `"abundance"`: Only predict cell-type proportions
  - `"expression"`: Predict both proportions and cell-type-specific expression
  - `"all"`: Predict proportions, expression, and perform DE analysis between conditions
  - `"all+"`: Same as "all" plus per-sample DE analysis (computationally intensive)
  - `"DE.only"`: Run DE analysis on existing CLIMB results
- **`up.lim`**: Upper bound for coefficients (L-infinity norm constraint). Default: `Inf`
- **`lambda`**: Regularization factor. Default: `0`
- **`verbose`**: Print progress messages. Default: `TRUE`
- **`cancer_pattern`**: String pattern to identify cancer cell-types for differential expression analysis. Default: `"*"`
- **`conditions`**: Vector of condition labels for DE analysis (required for "all" mode)
- **`min_common_genes`**: Minimum number of genes in common between bulk and single-cell datasets. Default: `100`
- **`ratio_cell_increase`**: Percentage increase at each empirical Bayes step. Default: `0.02`
- **`n.iter.subsampling`**: Number of subsampling iterations (results are averaged). Default: `5`
- **`min.n.cells`**: Minimum cells per cell-type to subsample. Default: `50`
- **`n.top_mean.genes`**: Number of genes for bulk-specific gene selection. Default: `500`

### Output

CLIMB returns a list containing:

- **`props.corrected`**: Matrix of cell-type proportions (samples × cell-types)
- **`expr.highres`**: Array of high-resolution cell-type expression (samples × genes × cell-types)
- **`expr.overall`**: Matrix of overall cell-type expression across all samples
- Additional elements depending on the mode (DE results, coefficients, etc.)

### Example

```r
# Load your data as ExpressionSet objects
library(Biobase)

# sc: single-cell reference with sc$cellType annotation
# bulk: bulk RNA-seq data to deconvolute

# Basic usage - predict proportions only
result <- climb(sc = sc, 
                bulk = bulk, 
                mode = "abundance")

# View proportions
head(result$props.corrected)

# Advanced usage - predict proportions and expression
result <- climb(sc = sc, 
                bulk = bulk, 
                mode = "expression",
                n.iter.subsampling = 10)

# View high-resolution expression for first sample, first 5 genes
result$expr.highres[1, 1:5, ]

# With differential expression analysis
result <- climb(sc = sc, 
                bulk = bulk, 
                mode = "all",
                conditions = bulk$condition,  # condition labels
                cancer_pattern = "Cancer")     # identify cancer cells
```

---

## CLIFF: Drug Sensitivity Deconvolution

### Description

CLIFF (Ce**L**l-type-specific drug sens**I**tivity **F**rom deconvoluted expression pro**F**iles) takes CLIMB output and integrates it with drug sensitivity data (AUC values) to identify cell-type-specific drug responses. It can optionally incorporate somatic mutation data for cancer cell-types.

### Key Features

- Estimates cell-type-specific drug sensitivity from bulk drug response data
- Integrates mutation data for cancer-specific effects
- Uses an expectation-maximization (EM) algorithm for robust estimation
- Supports both overall and high-resolution expression modes
- Optional L2 regularization (ridge regression)

### Usage

```r
cliff(climb_output, drug_data, mutation_data = NULL, 
      min.mutation = 0, max.em.steps = 100, 
      mode = "highres", regularization = "none", 
      cancer_pattern = "like")
```

### Parameters

- **`climb_output`**: Output object from CLIMB containing deconvoluted expression and proportions
- **`drug_data`**: Data frame with drug sensitivity data. Must contain:
  - `auc` column: Drug response AUC values
  - `sample` column: Sample identifiers matching CLIMB output
- **`mutation_data`**: Optional matrix (samples × mutations) with binary mutation status (1/0). Default: `NULL`
- **`min.mutation`**: Minimum number of samples with a mutation to include it in the model. Default: `0`
- **`max.em.steps`**: Maximum number of EM algorithm iterations. Default: `100`
- **`mode`**: Expression mode
  - `"highres"`: Use high-resolution deconvoluted expression from CLIMB (recommended)
  - `"overall"`: Use average cell-type expression across all samples
- **`regularization`**: Regularization method
  - `"none"`: No regularization (default)
  - `"L2"`: Apply ridge (L2) regularization
- **`cancer_pattern`**: String pattern to identify cancer cell-types for mutation attribution. Default: `"like"`

### Output

CLIFF returns a list with three elements:

1. **`PI_hat_nk`**: Matrix of cell-type-specific drug sensitivity predictions (samples × cell-types)
2. **`mutation_data`**: The mutation data matrix used in the analysis
3. **`climb_prop`**: Cell-type proportions from CLIMB

### Example

```r
# First run CLIMB to get deconvolution results
climb_result <- climb(sc = sc, 
                      bulk = bulk, 
                      mode = "expression")

# Prepare drug sensitivity data
drug_data <- data.frame(
  sample = colnames(bulk),
  auc = c(0.4, 0.6, 0.3, 0.7, 0.5, 0.8)  # AUC values
)

# Optional: prepare mutation data
mutation_data <- matrix(0, nrow = ncol(bulk), ncol = 5)
rownames(mutation_data) <- colnames(bulk)
colnames(mutation_data) <- c("TP53", "FLT3", "NPM1", "DNMT3A", "IDH1")
# Fill with 0/1 mutation status...

# Run CLIFF
cliff_result <- cliff(climb_output = climb_result,
                      drug_data = drug_data,
                      mutation_data = mutation_data,
                      mode = "highres",
                      cancer_pattern = "Cancer",
                      min.mutation = 3)

# View cell-type-specific drug sensitivity
head(cliff_result[[1]])  # PI_hat_nk matrix
```

---

## Workflow Example

Here's a complete workflow combining both methods:

```r
library(ClimbTheCliff)
library(Biobase)

# Step 1: Prepare your data
# sc: ExpressionSet with single-cell reference (must have sc$cellType)
# bulk: ExpressionSet with bulk RNA-seq data
# drug_data: data.frame with 'sample' and 'auc' columns
# mutation_data: matrix with mutation status (optional)

# Step 2: Run CLIMB deconvolution
climb_output <- climb(
  sc = sc,
  bulk = bulk,
  mode = "expression",
  verbose = TRUE,
  n.iter.subsampling = 10,
  min.n.cells = 50
)

# Step 3: Examine CLIMB results
print("Cell-type proportions:")
print(head(climb_output$props.corrected))

print("High-resolution expression dimensions:")
print(dim(climb_output$expr.highres))

# Step 4: Run CLIFF for drug sensitivity analysis
cliff_output <- cliff(
  climb_output = climb_output,
  drug_data = drug_data,
  mutation_data = mutation_data,
  mode = "highres",
  regularization = "L2",
  cancer_pattern = "Cancer"
)

# Step 5: Analyze cell-type-specific drug responses
cell_type_sensitivity <- cliff_output[[1]]
print("Cell-type-specific drug sensitivity:")
print(head(cell_type_sensitivity))

# Calculate weighted AUC per sample
weighted_auc <- rowSums(
  cell_type_sensitivity * climb_output$props.corrected
)
print("Predicted vs actual AUC:")
print(data.frame(
  actual = drug_data$auc,
  predicted = weighted_auc
))
```

---

## Citation

If you use ClimbTheCliff in your research, please cite:

*(Add citation information here once published)*

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Issues and Support

For bug reports, feature requests, or questions, please open an issue on the GitHub repository:

https://github.com/drayac/ClimbTheCliff/issues

---

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
