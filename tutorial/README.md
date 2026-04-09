# Reproducibility material for: *Optimal and computationally tractable lower bounds for logistic log-likelihoods*

This repository contains the code used to reproduce the numerical experiments presented in the article [Anceschi, et al. (2025)](https://arxiv.org/abs/2410.10309) and, more specifically, in Section 5 and Appendix F

In particular:
- Penalized likelihood experiments (Table 1 in the paper and supplement)
- Variational Bayes experiments (Figure 3 in the paper)

---

## Repository structure

Each subfolder corresponds to a specific application discussed in the paper and supplement.

```
tutorial/
├── Alzheimer/
├── Portland/
├── Spam/
└── tutorial_utils.R
```

The R script `tutorial_utils.R` contains shared utility functions used across all experiments, including plotting functions and the computation of accuracy metrics, such as the total variation distance between univariate distributions.


In particular, each application folder is organized as follows:

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

Each R script in the folder reproduces a specific part of the analysis discussed in the article for the corresponding application. Moreover, each R script is self-contained and ready-to-use. Specifically, they implement the follwing experiments:

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
2. Load the data `data/` folder from and standardize the predictors
3. Set all the optimization control parameters to default values
4. Run penalized likelihood optimization or variational inference under different bounds:
   - **BL** (Böhning–Lindsay)
   - **PG** (Polya-Gamma)
   - **PQ** (proposed method)
5. Save the outputs and results reported in the article or in the supplementary material:
   - `csv/` → numerical summaries (iterations, runtime, log-likelihood, etc.)
   - `img/` → optional diagnostic plots 
   - `rds/` → fitted models and intermediate results

The results reported in the article and supplementary material require the execution of all the R scripts under each of the application folders.
The following table provide a schematic representation of the main results in the paper and the scripts that can be used to reproduce them.

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

In particular, `Portland/portland_mfvb.R` typically takes a lot of time due to the expensive Gibb sampling procedure used to obtain a state-of-the-art estimate of the posterior distribution.

---

# Methodology 

- The implementation closely follows the methodology described in Section 5 of the paper and Sections D, E, and F of the Supplementary Material.
- All methods share the same optimization backbone, differing only in the choice of lower bound.
- The PQ method consistently reduces the number of iterations and often improves total runtime, especially in high-dimensional or non-smooth settings.

---

# References

- Niccolò Anceschi, Cristian Castiglione, Tommaso Rigon, Giacomo Zanella, Daniele Durante.  
  *Optimal and computationally tractable lower bounds for logistic log-likelihoods*  
  arXiv preprint: [arXiv:2410.10309](https://arxiv.org/abs/2410.10309)


