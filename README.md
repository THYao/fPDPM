# fPDPM: Functional Product Dirichlet Process Mixture

Research implementation of the functional product Dirichlet process mixture (fPDPM) model for two-dimensional functional/image observations, including wavelet-domain local clustering and a low-rank residual covariance component.

## Repository structure

| File | Description |
| --- | --- |
| `PDPM_Main.R` | Command-line entry point for fitting the model and saving posterior draws. |
| `Code/PDPM_dwt2d_LRD_func.R` | MCMC implementation, including the `fPDPM_mcmc()` function and supporting updates. |
| `Code/dmv_LRD_woodbury.cpp` | RcppArmadillo routine for a multivariate Gaussian log density with covariance `W W^T + s² I`, using the Woodbury identity. |
| `gwkConv.R` | Reconstructs the wavelet mean component for posterior draws and calculates location-wise Geweke diagnostics. |
| `postSummary_hclust.R` | Constructs posterior pairwise dissimilarities from local cluster memberships and performs hierarchical clustering. |

## Requirements

R packages used in the scripts include `MASS`, `extraDistr`, `MCMCpack`, `waveslim`, `Rcpp`, `RcppArmadillo`, `truncdist`, `matrixStats`, `stringr`, `ggplot2`, `coda`, and `cluster`.

Install missing dependencies in R, for example:

```r
install.packages(c(
  "MASS", "extraDistr", "MCMCpack", "waveslim", "Rcpp",
  "RcppArmadillo", "truncdist", "matrixStats", "stringr",
  "ggplot2", "coda", "cluster"
))
```

Compiling the C++ helper requires a working C++ toolchain and compatible linear-algebra/Fortran libraries on some platforms (including certain macOS R installations).

## Input data and directory setup

Create `Data/` and `MCMC_Res/` under the repository root. The input should be saved as `Data/<dataset>.RDS` and contain a list with an element named `df`:

```r
saveRDS(obsDf, "Data/example.RDS")
```

Here `obsDf` is a numeric array of dimensions **spatial row × spatial column × observation**. The current implementation assumes a square, power-of-two spatial grid for its 2D wavelet decomposition. 

**Important:** The three command-line R scripts contain `setwd("~/fPDPM/")`. Change that path to the local repository directory (or remove it and launch the scripts from the repository root) before running them.

## Simulated data

We provide the simulated dataset
[`diagErr_highSNR.RDS`](https://github.com/THYao/fPDPM/blob/main/Data/diagErr_highSNR.RDS)
to demonstrate the fPDPM implementation.

It contains **300 simulated observations on a 32 × 32 spatial grid**,
under a high signal-to-noise ratio setting with diagonal residual errors.

Load the dataset in R:

```r
obsDf <- readRDS("./Data/diagErr_highSNR.RDS")

dim(obsDf)
# [1] 32 32 300
```

The supplied scripts can be run with the dataset name
`diagErr_highSNR` and a chosen maximum factor
dimension (for example, `3`).



## Fit the model

The main MCMC function, `fPDPM_mcmc()`, is implemented in `Code/PDPM_dwt2d_LRD_func.R`. The script `PDPM_Main.R` provides an example workflow for fitting the model and saving the results.

### Running the example

We provide a simulated dataset,
[`diagErr_highSNR.RDS`](https://github.com/THYao/fPDPM/blob/main/Data/diagErr_highSNR.RDS),
containing 300 functional observations on a 32 × 32 spatial grid.

From the repository root directory, run:

```bash
Rscript PDPM_Main.R diagErr_highSNR 3 1
```

The two command-line arguments are:

1. `dfName`: Dataset filename without the `.RDS` extension.
2. `K_max`: Maximum number of factors in the low-rank covariance representation.

Alternatively, the model can be fitted directly in R:

```r
source("./Code/PDPM_dwt2d_LRD_func.R")

obsDf <- readRDS("./Data/diagErr_highSNR.RDS")

iteNum <- 3000
burnIn <- 2700
K_max <- 3
seed_ <- iteNum + K_max * 12

fit <- fPDPM_mcmc(
  obsDf = obsDf,
  iteNum = iteNum,
  burnIn = burnIn,
  wf_ = "haar",
  seed_ = seed_,
  K.max = K_max,
  KcutOff = 0.7,
  K.tol = 1e-1,
  saveCov = FALSE
)
```

The example runs 3,000 MCMC iterations, discards the first 2,700 as burn-in, and retains 300 posterior draws.

The implementation uses a two-dimensional discrete wavelet transform and a low-rank representation of the residual covariance. The number of factors is adaptively updated during MCMC.

**Note:** The provided `PDPM_Main.R` script uses the working directory `~/fPDPM/`. Update the `setwd()` statement to match your local repository location before running the script. Create an `MCMC_Res/` directory to store the results.

### MCMC output

The function `fPDPM_mcmc()` returns an R list containing:

| Component | Description |
|---|---|
| `initClu` | Initial cluster memberships |
| `llh` | Log-likelihood trace across all MCMC iterations |
| `clu` | Posterior cluster membership samples |
| `coef` | Posterior wavelet coefficient samples |
| `LRD` | Low-rank covariance parameters, if `saveCov = TRUE` |
| `var` | Posterior residual variance samples |

The main script saves the results as:

```text
MCMC_Res/<dfName>_Kmax<K_max>.RDS
```

The saved RDS object contains the MCMC execution time (`t`)
and fitted model results (`res`).

Posterior samples in `clu`, `coef`, and `var` are retained
only after burn-in. The log-likelihood trace (`llh`) includes
all iterations.


## Posterior mean component and Geweke diagnostics

```bash
Rscript gwkConv.R diagErr_highSNR 3 1
```

The script reconstructs the **wavelet mean component** from posterior cluster memberships and coefficients, selects every third retained draw, and computes Geweke Z-statistics for each observation and spatial location using `coda::geweke.diag()`. It converts each Z-statistic to a two-sided p-value using `2 * (1 - pnorm(abs(z)))`, applies Benjamini–Hochberg adjustment **within each observation across spatial locations**, and prints a summary of the proportion of adjusted p-values greater than 0.05.

**Current settings:** `gwkConv.R` hard-codes `nLoc = 32`, `nSmp = 300`, and `thinIn = 3`. Adapt these before analyzing other dimensions or sample sizes. 

## Pairwise dissimilarities and hierarchical clustering

```bash
Rscript postSummary_hclust.R diagErr_highSNR 3 1
```

For each retained posterior samples, the script computes pairwise dissimilarities based on agreement of local cluster memberships, excluding the final residual-covariance membership column. Each wavelet resolution has a specified weight. The resolution weights include the sample-size-dependent adjustment defined in the code and manuscript. Pairwise dissimilarities are averaged across draws, then passed to `hclust()`. The number of clusters is selected by maximizing average silhouette width over `k = 2, ..., 30`.

The script saves the per-draw dissimilarity array to `MCMC_Res/pairDist_<dataset>_Kmax<K>.RDS`, and plots the resulting dendrogram. Cluster assignments are calculated with `cutree()` but are not currently written to a file.

**Current settings:** `postSummary_hclust.R` hard-codes `nLoc = 32`, `nSmp = 300`, and `thinIn = 3`. Its `k` search up to 30 also assumes enough observations. Adapt these settings for other datasets.

## Notes on reproducibility

- The scripts are provided as research code, not as a packaged R library.
- Model-fitting runtime and MCMC mixing depend on the dataset and computational environment.
- Only post-burn-in samples are stored. 
- An example simulated dataset is provided in `Data/`. Users can also supply their own data following the specified input format.
- Precomputed MCMC results are not included. Create the `MCMC_Res/` directory before running the scripts.
- This repository has not been end-to-end tested in the current environment.

## Citation and license

If you use fPDPM in your research, please cite the accompanying methodological manuscript.

**Citation:** [Bayesian Nonparametric Global-local Clustering of High-dimensional Functions](https://arxiv.org/abs/2412.09792)

This repository is distributed under the MIT License.
See [LICENSE](LICENSE) for details.

