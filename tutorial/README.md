# Reproducibility material for ''Optimal and computationally tractable lower bounds for logistic log-likelihoods''

This repository contains the code required to reproduce the numerical experiments presented in the paper Anceschi, N., Castiglione, C., Rigon, T., Zanella, G. and Durante, D. (2025), [Optimal and computationally tractable lower bounds for logistic log-likelihoods](https://arxiv.org/abs/2410.10309), with particular focus on **Section 5** of the main paper and **Section F.2** of the Supplementary Material.

---

## Repository structure

The `tutorial/` folder is organized in sub-folders, each corresponding to a specific application in the paper and Supplementary Material:

- `Portland/` → Spatial modelling of motor-vehicle theft events in Portland (Oregon). [see, Section 5 and Supplementary Material F.2]
- `Alzheimer/` →Alzheimer’s disease classification based on demographic and biological covariates.  [see, Supplementary Material F.2]
- `Spam/` → Binary email classification (spam vs. non-spam) based on textual features. [see, Supplementary Material F.2]


Additionally, the `tutorial/` folder contains the R script `tutorial_utils.R`, which provides shared utility functions used across all experiments, including data preprocessing, plotting routines, and accuracy metrics (e.g., total variation distance).

Each application sub-folder is organized as follows:

- `csv/` → folder where numerical summaries (iterations, runtime, log-likelihood, etc.) will be saved
- `rds/` → folder where fitted models and intermediate results (in `.RDS` extension) will be saved
- `img/` → folder where optional diagnostic plots (Portland dataset only) will be saved
- `*_ridge.R` → Penalized estimation with **ridge** over a solution path [replace `*` with `portland`, `alzheimer` or `spam`]
- `*_lasso.R` → Penalized estimation with **lasso** over a solution path  [replace `*` with `portland`, `alzheimer` or `spam`]
- `*_enet.R` → Penalized estimation with **elastic-net** over a solution path  [replace `*` with `portland`, `alzheimer` or `spam`]
- `*_enet2d.R` → Penalized estimation with **elastic-net** over 2-dimensional tuning grid (only for Portland) [replace `*` with `portland`]
- `*_mfvb.R` → Bayesian inference via **variational Bayes** (only for Portland) [replace `*` with `portland`]

The above R scripts are self-contained and can be executed independently to reproduce specific parts of the analysis discussed in the article and the Supplementary Material.
All the application-specific experiments implemented in these R scripts follow the same structure:

1. Load the required packages and utility functions
2. Load the data from the `data/` folder and preprocess the predictors
3. Set the optimization control parameters (default values are used in all experiments)
4. Run penalized likelihood optimization or variational inference under different **tangent minorizers**:
   - **BL** (Böhning–Lindsay)
   - **PG** (Polya-Gamma)
   - **PQ** (piecewise quadratic, **proposed method**)
5. Save the diagnostic information and final results in the output folders `img/`, `csv/`, `rds/`

Reproducing the full set of results requires running all scripts across the three application folders (`Portland/`, `Alzheimer/`, and `Spam/`).
The following table provides a schematic representation of the main results in the paper and the scripts to reproduce them.

| Paper element | Location in repo |
|--------------|----------------|
| Table 1 (main paper) | `Portland/portland_enet.R`, `Portland/portland_enet2d.R` |
| Figure 2 (main paper) | `Portland/portland_mfvb.R` |
| Table F.1 (supplement) | all `*_ridge.R`, `*_lasso.R`, `*_enet.R` |

---

# How to run the codes and reproduce results

To run all experiments, do the following:
- set the working directory to the root folder of the package, i.e., the directory containing the `logitPQbound.Rproj` file.  Alternatively, RStudio users can simply open the `logitPQbound.Rproj` file, which will automatically set the working directory.  
- run separately each script [ `*_ridge.R`, `*_lasso.R`, `*_enet.R`,  `*_enet2d.R` (only for Portland), and `*_mfvb.R` (only for Portland), where `*` is replaced by `portland`, `alzheimer` or `spam` depending on the application considered]. For example, to run elastic-net under Portland data, write:

```r
# Set the working directory using the path to the package on your system
# setwd("~/logitPQbound/") 

source("tutorial/Portland/portland_enet.R")
```

As discussed in the article and Supplementary Material, the total runtime required for the different applications may range from a few minutes to several hours, depending on:

- the size of the dataset `n`
- the type of penalty `ridge`, `lasso`, `elastic-net`
- the number of regression parameters `p`
- the number of tuning parameters [1 for `ridge`, `lasso`, 2 for `elastic-net`]

---

# Reproducibility  notes

- The implementation closely follows the methodology described in **Section 5** of the paper and leverages the algorithms presented in detail **Sections D, E, and F** of the Supplementary Material.
- While the number of iterations to convergence is reproducible across systems, the total execution time depends on the specific hardware employed. Nonetheless, time variations apply to all the minorizers compared, and hence, the relative performance comparisons and final conclusions remain consistent across hardwares. The results reported in the paper and in the Supplementary Material were obtained on **MacBook Air (2020) equipped with an Apple M1 processor and 8 GB of RAM [R version 4.2.2]**.
---

# References

* Anceschi, N., Castiglione, C., Rigon, T., Zanella, G. and Durante, D. (2025)  
  [Optimal and computationally tractable lower bounds for logistic log-likelihoods](https://arxiv.org/abs/2410.10309), arXiv:2410.10309. 


