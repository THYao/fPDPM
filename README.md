# fPDPM: Functional Product Dirichlet Process Mixture

Research implementation of the functional product Dirichlet process mixture (fPDPM) model for two-dimensional functional/image observations, including wavelet-domain local clustering and a low-rank residual covariance component.

This repository includes a simulated dataset, model-fitting code, MCMC convergence diagnostics, and post-processing for global clustering.


## Repository structure

| File | Description |
| --- | --- |
| `PDPM_Main.R` | Command-line entry point for fitting the model and saving posterior draws. |
| `Code/PDPM_dwt2d_LRD_func.R` | MCMC implementation, including the `fPDPM_mcmc()` function and supporting updates. |
| `Code/dmv_LRD_woodbury.cpp` | RcppArmadillo routine for a multivariate Gaussian log density with covariance `W W^T + s² I`, using the Woodbury identity. |
| `chkConv.R` | Reconstructs observation-specific posterior quantities and calculates Geweke and Gelman–Rubin diagnostics. |
| `postSummary_hclust.R` | Computes weighted pairwise dissimilarities from local cluster assignments and performs hierarchical clustering. |
| `Data/diagErr_highSNR.RDS` | Example simulated observations (32 × 32 spatial locations for 300 observations). |
| `MCMC_Res/` | Local directory for MCMC results and clustering post-processing (create before running). |


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

Alternatively, the dataset can be downloaded directly from GitHub:

```r
url <- paste0("https://github.com/THYao/fPDPM/main/","Data/diagErr_highSNR.RDS")
download.file(
  url,
  destfile = "diagErr_highSNR.RDS",
  mode = "wb"
)

obsDf <- readRDS("diagErr_highSNR.RDS")
```

The supplied scripts can be run with the dataset name `diagErr_highSNR` and a chosen maximum factor dimension (for example, `3`).


## Fit the model

Create the results directory and run **two independent chains** from the repository root:

```bash
mkdir -p MCMC_Res
Rscript PDPM_Main.R diagErr_highSNR 3 1
Rscript PDPM_Main.R diagErr_highSNR 3 2
```

The three positional arguments are:

1. `dfName`: name of the input file, without `.RDS` (here `diagErr_highSNR`).
2. `K_max`: maximum factor dimension of the low-rank covariance representation (here `3`).
3. `chIdx`: chain index (here `1` or `2`). Each chain is assigned a different random seed.

The current `PDPM_Main.R` example runs **3,000** iterations per chain, discards the first **2,700** as burn-in, and saves **300** posterior draws. To change these settings, edit `iteNum` and `burnIn` in `PDPM_Main.R` before running. The longer 10,000-iteration real-data analysis described in the manuscript is a separate configuration, not the default in this example script.

Results are saved as:

```text
MCMC_Res/<dfName>_Kmax<K_max>_ch<chIdx>.RDS
```

For example:

```text
MCMC_Res/diagErr_highSNR_Kmax3_ch1.RDS
MCMC_Res/diagErr_highSNR_Kmax3_ch2.RDS
```

The main script calls `fPDPM_mcmc(..., saveCov = FALSE)`. The MCMC algorithm still updates its covariance parameters, but `saveCov = FALSE` omits their large low-rank (`LRD`) posterior objects from storage; **cluster memberships, wavelet coefficients, and variance components are still retained**.

## Saved MCMC output

Read one chain as follows:

```r
out <- readRDS("./MCMC_Res/diagErr_highSNR_Kmax3_ch1.RDS")
out$t       # Runtime information from system.time()
fit <- out$res
names(fit)
```

| Object | Contents |
| --- | --- |
| `fit$initClu` | Initial cluster membership matrix. |
| `fit$llh` | Log-likelihood trajectory across all iterations, including burn-in. |
| `fit$clu` | Array with dimensions observation × cluster-assignment block × retained iteration. The final block is for variance/covariance clustering; preceding blocks are for wavelet coefficients. |
| `fit$coef[[t]]` | At retained iteration `t`, a list of coefficient matrices, one per wavelet block. In each matrix, rows correspond to cluster indices, and columns contain the associated wavelet coefficients. |
| `fit$var[[t]]` | At retained iteration `t`, a vector of variance values indexed by the variance-cluster labels. |
| `fit$LRD[[t]]` | Low-rank covariance parameters, stored only when `saveCov = TRUE`. |

The objects `clu`, `coef`, and `var` contain **post-burn-in** draws only. The numerical cluster labels are iteration-specific and should not be interpreted as fixed observation identities across draws.

### Recovering the observation-specific mean and variance

The cluster-specific coefficients and variance parameters are **not themselves** the observation-specific trajectories used for convergence diagnostics. At each retained iteration `t`:

- For the wavelet mean component, use each observation's cluster assignments in the wavelet blocks to select the corresponding rows of `fit$coef[[t]]`; apply the inverse 2D wavelet transform to obtain the reconstructed function, denoted here by `theta_i(v)`.
- For the variance component, use the **last** column of the assignment matrix to select the variance assigned to each observation from `fit$var[[t]]`. Thus the observation-specific draw is `sigma2_i(t) = fit$var[[t]][s_i(t)]`.

These observation-specific quantities are invariant to permutations of the component labels when the assignments and component parameters are permuted together.

The model implementation already provides `gen_muMat()`, the same reconstruction used by `chkConv.R` (where the corresponding function is named `pdpm_muMat()`). This example retrieves **one retained posterior draw for every observation**:

```r
source("./Code/PDPM_dwt2d_LRD_func.R")
fit <- readRDS("./MCMC_Res/diagErr_highSNR_Kmax3_ch1.RDS")$res

t <- 1L                                 # First retained iteration
clu_t <- fit$clu[, , t]                 # n observations x membership blocks

theta_t <- gen_muMat(
  cluMat_ = clu_t,
  coefLs = fit$coef[[t]],
  wf_ = "haar"
)                                      # n observations x spatial locations

s_t <- clu_t[, ncol(clu_t)]            # Variance-cluster assignments
sigma2_t <- fit$var[[t]][s_t]          # One variance value per observation

# Values for observation i = 1:
theta_t[1, ]
sigma2_t[1]
```

To calculate posterior means of these observation-specific quantities without storing all reconstructed functions at once:

```r
n_draw <- dim(fit$clu)[3]
n_obs <- dim(fit$clu)[1]
n_spatial <- ncol(theta_t)

theta_sum <- matrix(0, nrow = n_obs, ncol = n_spatial)
sigma2_sum <- numeric(n_obs)

for (t in seq_len(n_draw)) {
  clu_t <- fit$clu[, , t]
  theta_sum <- theta_sum + gen_muMat(clu_t, fit$coef[[t]], wf_ = "haar")
  sigma2_sum <- sigma2_sum + fit$var[[t]][clu_t[, ncol(clu_t)]]
}

theta_hat <- theta_sum / n_draw      # Posterior mean of reconstructed functions
sigma2_hat <- sigma2_sum / n_draw    # Posterior mean of assigned variances
```

For **trace plots and convergence diagnostics**, use the corresponding value at **every retained iteration**, not just these posterior averages. The next section shows how `chkConv.R` constructs and checks those trajectories.

## MCMC convergence diagnostics (`chkConv.R`)

First fit and save **both chains** using the commands above. Then run:

```bash
Rscript chkConv.R diagErr_highSNR 3 1
```

The arguments identify the dataset, `K_max`, and the chain index used for the **Geweke** calculations. The script's **Gelman–Rubin** section currently sets `nChain <- 2` and reads chain files `_ch1.RDS` and `_ch2.RDS`, irrespective of the third argument. Both files must exist.

**Current example-specific settings:** `chkConv.R` sets `nLoc = 32`, `nSmp = 300`, and `thinIn = 3`. Change these values for other datasets or sampling intervals. The thinning code selects the 3rd, 6th, 9th, ... retained draws; for 300 retained draws, the diagnostics use 100 per chain.

### What is reconstructed

At each selected retained iteration, `chkConv.R` builds:

| Object | Dimensions | Meaning |
| --- | --- | --- |
| `muArr` | observations × spatial locations × selected draws | Reconstructed wavelet mean, `theta_i(v)`, for the chain used in Geweke diagnostics. |
| `s2Mat` | observations × selected draws | Observation-specific assigned variance, `sigma2_{s_i}`, for that chain. |
| `estLs[[chIdx]]$mu` | observations × spatial locations × selected draws | Reconstructed mean trajectories for chain `chIdx`, used for Gelman–Rubin. |
| `estLs[[chIdx]]$s2` | observations × selected draws | Assigned-variance trajectories for chain `chIdx`, used for Gelman–Rubin. |

The variance draws are retrieved with the observation's **final-column membership**, rather than treating an arbitrarily numbered DP variance component as an identifiable fixed parameter.

### Geweke diagnostics (within one chain)

The script computes:

- `gwkMu.arr[i, v]`: Geweke Z-statistic for the reconstructed mean at observation `i`, spatial location `v`.
- `gwkS2.vec[i]`: Geweke Z-statistic for the assigned variance of observation `i`.
- `fdr_p_Prop[i]`: proportion of spatial locations for observation `i` whose **BH-adjusted two-sided Geweke p-values** exceed 0.05.

The current script prints `summary(fdr_p_Prop)` and `mean(gwkS2.vec)`. **The latter is the mean Geweke Z-statistic, not a p-value or a convergence proportion.** To inspect variance p-values, add this to the end of the Geweke section or run it in the same R session after the objects are constructed:

```r
gwkS2.p <- 2 * pnorm(-abs(gwkS2.vec))
summary(gwkS2.p)
mean(gwkS2.p > 0.05, na.rm = TRUE)
```

### Gelman–Rubin diagnostics (two chains)

The script uses `coda::gelman.diag(..., autoburnin = FALSE, multivariate = FALSE)` and compares the **same observations and spatial locations** across the two independently fitted chains.

| Object | Dimensions | Meaning |
| --- | --- | --- |
| `GR_mu[i, v]` | observations × spatial locations | Gelman–Rubin point estimate for the reconstructed wavelet mean at observation `i`, location `v`. |
| `GR_s2[i, 1]` | observations × 1 | Gelman–Rubin point estimate for observation-specific assigned variance. |
| `GR_s2[i, 2]` | observations × 1 | Upper confidence limit reported by `coda::gelman.diag()` for the corresponding variance diagnostic. |

The script currently prints summaries and proportions below the threshold 1.2:

```r
summary(as.vector(GR_mu))
mean(GR_mu < 1.2, na.rm = TRUE)
summary(GR_s2[, 1])
mean(GR_s2[, 1] < 1.2, na.rm = TRUE)
```

These diagnostics describe agreement and mixing for the **reconstructed observation-specific quantities**, not every latent component of the full covariance model. Interpret undefined or non-finite diagnostics separately rather than counting them as converged.

**Saving diagnostic outputs (optional):** The current `chkConv.R` **prints** these summaries but does not automatically write the diagnostic matrices to a file. If you want to retain them, append this block at the end of `chkConv.R`:

```r
saveRDS(
  list(
    geweke_mu_z = gwkMu.arr,
    geweke_sigma2_z = gwkS2.vec,
    geweke_mu_BH_nonsignificant_proportion = fdr_p_Prop,
    GR_mu = GR_mu,
    GR_sigma2 = GR_s2
  ),
  file = paste0("./MCMC_Res/", dfName, "_Kmax", K_max,
                "_ch", as.integer(args[3]), "_diagnostics.RDS")
)
```

## Global clustering (`postSummary_hclust.R`)

Run the post-processing script after fitting a chain:

```bash
Rscript postSummary_hclust.R diagErr_highSNR 3 1
```

The script uses the **wavelet-block** cluster assignments, omitting the final column used for variance/covariance clustering. For each selected posterior draw it calculates weighted pairwise dissimilarities between observations based on agreement of local cluster memberships. It averages these dissimilarities across draws, applies hierarchical clustering with `hclust()`, and selects the number of global clusters by maximizing average silhouette width over `k = 2, ..., 30`.

The per-draw pairwise dissimilarity array is saved as:

```text
MCMC_Res/pairDist_<dfName>_Kmax<K_max>_ch<chIdx>.RDS
```

The final cluster assignments are computed by `cutree()`, but the current script does not save them. To reproduce those assignments from the saved dissimilarity array:

```r
library(cluster)

pairDist <- readRDS(
  "./MCMC_Res/pairDist_diagErr_highSNR_Kmax3_ch1.RDS"
)
distMat <- apply(pairDist, c(1, 2), mean)
hc <- hclust(as.dist(distMat))

k_grid <- 2:30
sil_score <- vapply(k_grid, function(k) {
  summary(silhouette(cutree(hc, k = k), dist = as.dist(distMat)))$avg.width
}, numeric(1))

k_opt <- k_grid[which.max(sil_score)]
global_clusters <- cutree(hc, k = k_opt)

k_opt
global_clusters
plot(hc)
```

**Current example-specific settings:** `postSummary_hclust.R` uses `nLoc = 32`, `nSmp = 300`, and `thinIn = 3`. Adjust these settings as needed; the silhouette search over `2:30` assumes at least 30 observations.

## Notes on reproducibility

- The scripts are provided as research code, not as a packaged R library.
- The example dataset is included under `Data/`; users can also provide their own array in the required format.
- Precomputed MCMC chain outputs are not included. Create `MCMC_Res/` and run the model before post-processing.
- Fitting time and MCMC behavior depend on the dataset, hyperparameters, number of iterations, and computational environment.
- `saveCov = FALSE` reduces storage by omitting low-rank covariance draws while retaining the values needed to reconstruct the wavelet mean and observation-specific variance trajectories.

## Citation and license

If you use fPDPM in your research, please cite the accompanying methodological manuscript.

**Citation:** [Bayesian Nonparametric Global-local Clustering of High-dimensional Functions](https://arxiv.org/abs/2412.09792)

This repository is distributed under the MIT License.
See [LICENSE](LICENSE) for details.

