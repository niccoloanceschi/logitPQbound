#' @keywords internal
"_PACKAGE"

#' @useDynLib logitPQbound, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom Matrix crossprod
#' @importFrom Matrix Cholesky
#' @importFrom Matrix solve
#' @importFrom Matrix t
NULL