
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <chrono>

using namespace Rcpp;
using namespace arma;

//'  accurate inverse-matrix vector product Q^{-1}r via Cholesky factorization
//'
//' @param cholQ the Cholesky factor of Q
//' @param r A vector
//'
//' @export
// [[Rcpp::export]]
arma::vec invQ_r_Rcpp(const arma::mat& cholQ, const arma::vec& r) {
  return solve(trimatu(cholQ), solve(trimatl(cholQ.t()), r));
}

//' projected gradient descent, optimized for large-p-small-n problems
//'
//' @param u (arma::vec(n)) initial solution
//' @param y_nu (arma::vec(n)) auxiliary vector 1
//' @param nu_w (arma::vec(n)) auxiliary vector 2
//' @param xPxt_nu (arma::mat(n,n)) auxiliary
//' @param cholQ (arma::mat(n,n)) auxiliary Cholesky matrix
//' @param alpha (double) gradient step size
//' @param thrs (double) stopping threshold on relative absolute change
//' @param maxiter (int) maximum number of iterations
//' 
//' @return logistic penalized log likelihood
//' 
//' @export
// [[Rcpp::export]]
arma::vec prjgrad_large_p(arma::vec u, 
                          const arma::vec& y_nu, 
                          const arma::vec& nu_w, 
                          const arma::mat& xPxt_nu, 
                          const arma::mat& cholQ, 
                          double alpha, 
                          double thrs,
                          int maxiter) {
  arma::vec uOld = u;
  double u_diff = 2*thrs;
  int l = 0;
  while (l < maxiter && u_diff > thrs) {
    u -= alpha * (xPxt_nu.t() * invQ_r_Rcpp(cholQ, nu_w % (u - y_nu)));
    u = clamp(u, -1.0, 1.0);
    u_diff = arma::max(arma::abs(u - uOld) / (arma::abs(uOld) + 1e-8));
    uOld = u;
    l += 1;
  }
  return u;
}

//' projected gradient descent, optimized for large-n-small-p problems
//'
//' @param u (arma::vec(n)) initial solution
//' @param y_nu (arma::vec(n)) auxiliary vector
//' @param nu_x (arma::mat(n,p)) auxiliary matrix
//' @param cholQ (arma::mat(p,p)) auxiliary Cholesky matrix
//' @param alpha (double) gradient step size
//' @param thrs (double) stopping threshold on relative absolute change
//' @param maxiter (int) maximum number of iterations
//' 
//' @return logistic penalized log likelihood
//' 
//' @export
// [[Rcpp::export]]
arma::vec prjgrad_large_n(arma::vec u, 
                          const arma::vec& y_nu, 
                          const arma::mat& nu_x, 
                          const arma::mat& cholQ, 
                          double alpha, 
                          double thrs,
                          int maxiter) {
  arma::vec uOld = u;
  double u_diff = 2*thrs;
  int l = 0;
  while (l < maxiter && u_diff > thrs) {
    u -= alpha * (nu_x * invQ_r_Rcpp(cholQ, nu_x.t() * (u - y_nu)));
    u = clamp(u, -1.0, 1.0); 
    u_diff = arma::max(arma::abs(u - uOld) / (arma::abs(uOld) + 1e-8));
    uOld = u;
    l += 1;
  }
  return u;
}


//' inverse logit function
//'
//' @param eta (arma::vec(n)) linear predictors
//' 
//' @return inverse logit
//' 
//' @export
// [[Rcpp::export]]
arma::vec plogis(arma::vec eta) {
  return(1/(1+exp(-eta)));
}

//' logistic penalized log likelihood under elastic net regularization
//'
//' @param eta (arma::vec(n)) linear predictors
//' @param y (arma::vec(n)) responses
//' @param beta (arma::vec(p)) parameters (including intercept)
//' @param alpha (double) relative weight of L1 and L2 penalty
//' @param lambda (double) regularization parameter
//' 
//' @return logistic penalized log likelihood
//' 
//' @export
// [[Rcpp::export]]
double logpost(arma::vec eta, arma::vec y, arma::vec beta, double alpha, double lambda) {
  int p = beta.n_elem;
  arma::vec beta_pen  = beta.subvec(1,p-1);
  return sum(y % eta - log(1 + exp(eta))) - lambda*(alpha*sum(abs(beta_pen)) + (1-alpha)*0.5*sum(pow(beta_pen,2)));
}

//' PG bound coefficients
//'
//' @param eta (arma::vec(n)) linear predictors
//' 
//' @return PG weights
//' 
//' @export
// [[Rcpp::export]]
arma::vec w_PG(arma::vec eta) {
  arma::vec out(eta.n_elem);
  out = tanh(eta/2)/(2*eta);
  out.replace(arma::datum::nan,0.25);
  return(out);
}


//' PQ bound coefficients
//'     Remark: correcting for numerical instability around 0
//'             numerical instability still present at Infinity  
//'
//' @param eta (arma::vec(n)) linear predictors
//' 
//' @return PQ w and nu (as an arma::vec(2*n))
//'
//' @export
// [[Rcpp::export]]
arma::vec coeff_PQ(arma::vec eta) {

  arma::vec etaA = arma::abs(eta);
  arma::vec eta2 = arma::pow(eta,2);

  arma::vec vTanh = arma::tanh(eta/2)/(2*eta);
  arma::vec vLogCosh = 2*arma::log(arma::cosh(eta/2))/eta2; // [version 1]
  // arma::vec vLogCosh = 1/etaA + 2*(arma::log1p(arma::exp(-etaA))-log(2.))/eta2; // [version 2]

  arma::vec w_out = -vLogCosh+2.*vTanh;
  arma::vec nu_out = etaA%(vLogCosh-vTanh);

  arma::uvec eta0 = arma::find(etaA<1.0e-3); 
  
  w_out(eta0) = 0.25 - eta2(eta0)/32.;
  nu_out(eta0) = etaA(eta0) % eta2(eta0)/96.;

  return(arma::join_cols(w_out,nu_out));
}


