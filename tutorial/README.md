# Reproducibility material for: *Optimal and computationally tractable lower bounds for logistic log-likelihoods*

This repository contains the code required to reproduce the numerical experiments presented in [Anceschi, et al. (2025)](https://arxiv.org/abs/2410.10309), with particular focus on Section 5 of the main paper and Section F of the supplementary material.

---

## Repository structure

The repository is organized in sub-folders, each corresponding to a specific application discussed in the paper and supplement:

```
tutorial/
├── Alzheimer/
├── Portland/
├── Spam/
└── tutorial_utils.R
```

The R script `tutorial_utils.R` provides shared utility functions used across all experiments, including data preprocessing, plotting routines, and accuracy metrics (e.g., total variation distance).

Each application folder is organized as follows:

```
tutorial/
├── <Application>/
│ ├── csv/
│ ├── img/
│ ├── rds/
│ ├── <application>_enet.R
│ ├── <application>_lasso.R
│ ├── <application>_ridge.R
│ └── <application>_mfvb.R
├── ...
```

Inside each application folder, all the R scripts are self-contained and can be executed independently to reproduce specific parts of the analysis discussed in the article.
Specifically, each script implements one of the following experiments settings:

- `<application>_ridge.R`  
  → Penalized likelihood with **ridge penalty**  
- `<application>_lasso.R`  
  → Penalized likelihood with **lasso penalty**  
- `<application>_enet.R`  
  → Penalized likelihood with **elastic-net penalty**
- `<application>_mfvb.R`  
  → Bayesian inference via **variational Bayes** (only for Portland data)

All the experiments listed above follow the same structure:

1. Load the required packages and the utility functions
2. Load the data from the `data/` folder and standardize the predictors
3. Set the optimization control parameters (default values are used in all experiments)
4. Run penalized likelihood optimization or variational inference under different bounds:
   - **BL** (Böhning–Lindsay)
   - **PG** (Polya-Gamma)
   - **PQ** (proposed method)
5. Save the outputs required to reproduce the results reported in the paper and supplementary material:
   - `csv/` → numerical summaries (iterations, runtime, log-likelihood, etc.)
   - `img/` → optional diagnostic plots 
   - `rds/` → fitted models and intermediate results

Reproducing the full set of results requires running all scripts across the three application folders.
The following table provides a schematic representation of the main results in the paper and the scripts that can be used to reproduce them.

| Paper element | Location in repo |
|--------------|----------------|
| Table 1 (main paper) | `Portland/portland_enet.R` |
| Figure 3 (main paper) | `Portland/portland_mfvb.R` |
| Table F.1 (supplement) | all `*_ridge.R`, `*_lasso.R`, `*_enet.R` |

---

# How to reproduce results

R users can easily run each script independently using the following example code:

```r
source("tutorial/Portland/portland_enet.R")
```

As detailed in the article and supplementary material, the numerical experiments may take from minutes to hours depending on:

- dataset size
- penalty type
- number of regression parameters
- number of tuning parameters

In particular, `Portland/portland_mfvb.R` typically takes a lot of time due to the expensive Gibbs sampling procedure used to obtain a state-of-the-art estimate of the posterior distribution.

---

# Methodology 

- The implementation closely follows the methodology described in Section 5 of the paper and Sections D, E, and F of the supplementary material.
- The PQ method consistently reduces the number of iterations and often improves total runtime, especially in high-dimensional or non-smooth settings.

---

# References

- Niccolò Anceschi, Cristian Castiglione, Tommaso Rigon, Giacomo Zanella, Daniele Durante.  
  *Optimal and computationally tractable lower bounds for logistic log-likelihoods*  
  arXiv preprint: [arXiv:2410.10309](https://arxiv.org/abs/2410.10309)


