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

- `img/` → optional diagnostic plots
- `csv/` → numerical summaries (iterations, runtime, log-likelihood, etc.)
- `rds/` → fitted models and intermediate results (in `.RDS` extension)
- `*_logit.R` → Un-penalized likelihood (Spam dataset only)
- `*_ridge.R` → Penalized likelihood with **ridge penalty** over a solution path
- `*_lasso.R` → Penalized likelihood with **lasso penalty** over a solution path
- `*_enet.R` → Penalized likelihood with **elastic-net penalty** over a solution path
- `*_enet2d.R` → Penalized likelihood with **elastic-net penalty** over a two-dimensional tuning grid (Portland dataset only)
- `*_mfvb.R` → Bayesian inference via **variational Bayes** (Portland dataset only)

The above R scripts are self-contained and can be executed independently to reproduce specific parts of the analysis discussed in the article.
All the application-specific experiments implemented in such R scripts follow the same structure:

1. Load the required packages and utility functions
2. Load the data from the `data/` folder and preprocess the predictors
3. Set the optimization control parameters (default values are used in all experiments)
4. Run penalized likelihood optimization or variational inference under different MM bounds:
   - **BL** (Böhning–Lindsay)
   - **PG** (Polya-Gamma)
   - **PQ** (piecewise quadratic, proposed method)
5. Save the diagnostic information and final results in the output folders `img/`, `csv/`, `rds/`

Reproducing the full set of results requires running all scripts across the three application folders (`Portland/`, `Alzheimer/`, and `Spam/`).
The following table provides a schematic representation of the main results in the paper and the scripts that can be used to reproduce them.

| Paper element | Location in repo |
|--------------|----------------|
| Table 1 (main paper) | `Portland/portland_enet.R`, `Portland/portland_enet2d.R` |
| Figure 3 (main paper) | `Portland/portland_mfvb.R` |
| Table F.1 (supplement) | all `*_ridge.R`, `*_lasso.R`, `*_enet.R` |

---

# How to reproduce results

To run all experiments, it is necessary to set the working directory to the root folder of the package, i.e., the directory containing the `logitPQbound.Rproj` file.  
Alternatively, RStudio users can simply open the `logitPQbound.Rproj` file, which will automatically set the working directory.  

R users can run each script independently using, for example:

```r
# Set the working directory using the path to the package on your system
# setwd("~/logitPQbound/") 

source("tutorial/Portland/portland_enet.R")
```

As discussed in the article and supplementary material, the computational time required for the numerical experiments may range from a few minutes to several hours, depending on:

- the size of the dataset
- the type of penalty
- the number of regression parameters
- the number of tuning parameters

---

# Methodology 

- The implementation closely follows the methodology described in Section 5 of the paper and Sections D, E, and F of the supplementary material.
- The PQ method consistently reduces the number of iterations and almost always improves total runtime, especially in high-dimensional or non-smooth settings.
- While the number of iterations to convergence is reproducible across systems, the execution time depends on the specific hardware. In particular, the results reported in the paper were obtained on a Dell XPS 15 laptop equipped with a 4.7 GHz processor and 32 GB of RAM.

---

# References

* Anceschi, N. Castiglione, C. and Rigon, T. and Zanella, G. and Durante, D. (2025)  
  [Optimal and computationally tractable lower bounds for logistic log-likelihoods](https://arxiv.org/abs/2410.10309), ArXiv. 


