# logitPQbound

The repository includes the R functions to implement the novel piece-wise qaudratic bound for logistic loglikelihoods from the paper

* Anceschi, N. and Rigon, T. and Zanella, G. and Durante, D. (2024), [Optimal lower bounds for logistic log-likelihoods](https://arxiv.org/abs/2410.10309), ArXiv. 

The `logitPQbound` package can be installed by running the following commands

```R
# If the devtools R package is not already installed
# install.packages("devtools")

devtools::install_github("niccoloanceschi/logitPQbound")
```
## Structure

- `src/` contains the Rccp code scripts
- `R/` contains the R source scripts
- `data/` contains input data
- `turorial/` contains code to reproduce results and plots from main paper

## Requirementss

Install packages:
```R
install.packages(c("devtools", "RcppArmadillo", "ggplot2", "dplyr"))
```
