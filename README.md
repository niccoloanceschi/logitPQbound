
# logitPQbound

The repository provides the R functions implementing the novel piecewise quadratic bound for logistic log-likelihoods from the paper

* Anceschi, N., Castiglione, C., Rigon, T., Zanella, G., and Durante, D. (2025), [Optimal and computationally tractable lower bounds for logistic log-likelihoods](https://arxiv.org/abs/2410.10309), arXiv:2410.10309. 

The `logitPQbound` package can be installed by running the following R commands

```R
# If the devtools R package is not already installed
# install.packages("devtools")

devtools::install_github("niccoloanceschi/logitPQbound")
```

Alternatively, if `devtools::install_github()` produces warnings or fails, the `pak` package can be used:

```R
# If the pak R package is not already installed
# install.packages("pak")

pak::pak("niccoloanceschi/logitPQbound")
```

## Structure

- `src/` contains the Rccp code scripts
- `R/` contains the R source scripts
- `man/` contains the documentation files for the exported R functions (generated via `roxygen2`)
- `data/` contains example datasets included in the package (in `.RData` extension)
- `turorial/` contains scripts to reproduce the numerical results and figures from the main paper
- `DESCRIPTION` is the package metadata file specifying authors, maintainers, version, and dependencies 
- `NAMESPACE` defines both the imported and exported functions
- `LICENSE` contains the licensing terms of the package (MIT license)  
- `.Rbuildignore` specifies files and directories that should be excluded when building the package  
- `.gitignore` specifies files and directories ignored by Git version control  
- `logitPQbound.Rproj` is the RStudio project file for development and reproducibility  

## Requirements

Install packages:

```R
install.packages(c("devtools", "pak", "RcppArmadillo", "Rmpfr", "Matrix", "sparseinv", "BayesLogit", "ggplot2", "dplyr"))
```
