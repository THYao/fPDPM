# fPDPM: Functional Product Dirichlet Process Mixture

Research implementation of the functional product Dirichlet process mixture (fPDPM) model for two-dimensional functional/image observations, including wavelet-domain local clustering and a low-rank residual covariance component.

## Repository structure

| File | Description |
| --- | --- |
| `PDPM_Main.R` | Command-line entry point for fitting the model and saving posterior draws. |
| `Code/PDPM_dwt2d_LRD_func.R` | MCMC implementation, including the `fPDPM_mcmc()` function and supporting updates. |
| `Code/dmv_LRD_woodbury.cpp` | RcppArmadillo routine for a multivariate Gaussian log density with covariance `W W^T + s² I`, using the Woodbury identity. |
| `Code/gwkConv.R` | Reconstructs the wavelet mean component for posterior draws and calculates location-wise Geweke diagnostics. |
| `Code/postSummary_hclust.R` | Constructs posterior pairwise dissimilarities from local cluster memberships and performs hierarchical clustering. |

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
saveRDS(list(df = obsDf), "Data/example.RDS")
```

Here `obsDf` is a numeric array of dimensions **spatial row × spatial column × observation**. The current implementation assumes a square, power-of-two spatial grid for its 2D wavelet decomposition. **No example data are supplied in this repository.**

**Important:** The three command-line R scripts contain `setwd("~/fPDPM/")`. Change that path to the local repository directory (or remove it and launch the scripts from the repository root) before running them.

## Simulated data

We provide the simulated dataset
[`diagErr_highSNR_nLoc32_nSmp300_b3.RDS`](https://github.com/THYao/fPDPM/blob/main/Data/diagErr_highSNR_nLoc32_nSmp300_b3.RDS)
to demonstrate the fPDPM implementation.

It contains **300 simulated observations on a 32 × 32 spatial grid**,
under a high signal-to-noise ratio setting with diagonal residual errors.

Load the dataset in R:

```r
obsDf <- readRDS("./Data/diagErr_highSNR_nLoc32_nSmp300_b3.RDS")$df

dim(obsDf)
# [1] 32 32 300
```

The supplied scripts can be run with the dataset name
`diagErr_highSNR_nLoc32_nSmp300_b3` and a chosen maximum factor
dimension (for example, `3`).


## Fit the model

```bash
Rscript PDPM_Main.R example 3
```

The first argument is the dataset filename without `.RDS`; the second is `K.max`, the maximum number of low-rank factors. The supplied script uses 3,000 iterations, 2,700 burn-in iterations, the Haar wavelet, `KcutOff = 0.7`, `K.tol = 0.1`, and `saveCov = FALSE`. It saves:

```text
MCMC_Res/example_Kmax3.RDS
```

The saved object has elements `t` (elapsed time) and `res` (MCMC results). The main model interface is:

```r
fPDPM_mcmc(
  obsDf, iteNum, burnIn, seed_ = 1, wf_,
  K.max = NULL, KcutOff = NULL, K.tol = 1e-1,
  init.nClu = 3, hyperAlpha = 3, saveCov = TRUE
)
```

The returned list includes `initClu`, `llh`, `clu`, `coef`, `LRD`, and `var`. `clu` has dimensions **observation × membership component × retained iteration**. The final membership column corresponds to the residual covariance/variance component; the other columns correspond to the three detail subbands (LH, HL, HH) at each wavelet resolution. When `saveCov = FALSE`, low-rank covariance draws are not retained.

## Posterior mean component and Geweke diagnostics

```bash
Rscript gwkConv.R example 3
```

The script reconstructs the **wavelet mean component** from posterior cluster memberships and coefficients, selects every third retained draw, and computes Geweke Z-statistics for each observation and spatial location using `coda::geweke.diag()`. It converts each Z-statistic to a two-sided p-value using `2 * (1 - pnorm(abs(z)))`, applies Benjamini–Hochberg adjustment **within each observation across spatial locations**, and prints a summary of the proportion of adjusted p-values greater than 0.05.

**Current settings:** `gwkConv.R` hard-codes `nLoc = 32`, `nSmp = 300`, and `thinIn = 3`. Adapt these before analyzing other dimensions or sample sizes. 

## Pairwise dissimilarities and hierarchical clustering

```bash
Rscript postSummary_hclust.R example 3
```

For each retained posterior samples, the script computes pairwise dissimilarities based on agreement of local cluster memberships, excluding the final residual-covariance membership column. Each wavelet resolution has a specified weight. The resolution weights include the sample-size-dependent adjustment defined in the code and manuscript. Pairwise dissimilarities are averaged across draws, then passed to `hclust()`. The number of clusters is selected by maximizing average silhouette width over `k = 2, ..., 30`.

The script saves the per-draw dissimilarity array to `MCMC_Res/pairDist_<dataset>_Kmax<K>.RDS`, and plots the resulting dendrogram. Cluster assignments are calculated with `cutree()` but are not currently written to a file.

**Current settings:** `postSummary_hclust.R` hard-codes `nLoc = 32`, `nSmp = 300`, and `thinIn = 3`. Its `k` search up to 30 also assumes enough observations. Adapt these settings for other datasets.

## Notes on reproducibility

- The scripts are provided as research code, not as a packaged R library.
- Model-fitting runtime and MCMC mixing depend on the dataset and computational environment.
- Only post-burn-in samples are stored. 
- No dataset or precomputed results are included; create the required directories and supply your own data.
- This repository has not been end-to-end tested in the current environment.

## Citation and license

Please add the manuscript citation and choose a license before public release. No license is granted by default.
