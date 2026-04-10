# Reproducibility material for: *Optimal and computationally tractable lower bounds for logistic log-likelihoods*

This repository contains the code required to reproduce the numerical experiments presented in the paper Anceschi, N. and Castiglione, C. and Rigon, T. and Zanella, G. and Durante, D. (2025), [Optimal and computationally tractable lower bounds for logistic log-likelihoods](https://arxiv.org/abs/2410.10309), with particular focus on Section 5 of the main paper and Section F of the supplementary material.

---

## Repository structure

The `tutorial/` folder is organized in sub-folders, each corresponding to a specific application discussed in the paper and supplementary material:

- `Portland/` → Spatial risk modelling of motor-vehicle theft events in Portland (Oregon).
- `Alzheimer/` → Early-stage Alzheimer’s disease classification based on demographic and biological covariates.
- `Spam/` → Binary email classification (spam vs. non-spam) based on textual features.

Additionally, the `tutorial/` folder contains the R script `tutorial_utils.R`, which provides shared utility functions used across all experiments, including data preprocessing, plotting routines, and accuracy metrics (e.g., total variation distance).

Each application sub-folder is organized as follows:

- `img/`
- `csv/`
- `rds/`
- `*_enet.R`
- `*_lasso.R`
- `*_ridge.R`
- `*_mfvb.R`

The above R scripts are self-contained and can be executed independently to reproduce specific parts of the analysis discussed in the article.
Specifically, each script implements one of the following experiments settings:

- `*_ridge.R` → Penalized likelihood with **ridge penalty**  
- `*_lasso.R` → Penalized likelihood with **lasso penalty**  
- `*_enet.R` → Penalized likelihood with **elastic-net penalty**
- `*_mfvb.R` → Bayesian inference via **variational Bayes** (only for Portland data)

All these experiments follow the same structure:

1. Load the required packages and utility functions
2. Load the data from the `data/` folder and preprocess the predictors
3. Set the optimization control parameters (default values are used in all experiments)
4. Run penalized likelihood optimization or variational inference under different MM bounds:
   - **BL** (Böhning–Lindsay)
   - **PG** (Polya-Gamma)
   - **PQ** (piecewise quadratic, proposed method)
5. Save the outputs and the results reported in the paper and supplementary material:
   - `img/` → optional diagnostic plots 
   - `csv/` → numerical summaries (iterations, runtime, log-likelihood, etc.)
   - `rds/` → fitted models and intermediate results

Reproducing the full set of results requires running all scripts across the three application folders (`Portland/`, `Alzheimer/`, and `Spam/`).
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

---

# Methodology 

- The implementation closely follows the methodology described in Section 5 of the paper and Sections D, E, and F of the supplementary material.
- The PQ method consistently reduces the number of iterations and often improves total runtime, especially in high-dimensional or non-smooth settings.

---

# References

* Anceschi, N. Castiglione, C. and Rigon, T. and Zanella, G. and Durante, D. (2025), [Optimal and computationally tractable lower bounds for logistic log-likelihoods](https://arxiv.org/abs/2410.10309), ArXiv. 


