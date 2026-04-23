
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <chrono>

using namespace Rcpp;
using namespace arma;

//' PG bound coefficients
//'
//' @param eta (arma::vec(n)) linear predictors
//' 
//' @return PG weights
//' 
//' @export
// [[Rcpp::export]]
arma::vec coeff_PG(arma::vec eta) {
  arma::vec out(eta.n_elem);
  out = arma::tanh(eta/2)/(2*eta);
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
arma::mat coeff_PQ(arma::vec eta) {
  
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
  
  return(arma::join_rows(w_out,nu_out));
}
 
arma::vec chol_solve(const arma::mat& cholQ, const arma::vec& r) {
  return solve(trimatu(cholQ), solve(trimatl(cholQ.t()), r));
}

arma::vec soft_threshold(const arma::vec& x, const double& lambda) {
  arma::vec z = arma::abs(x) - lambda;
  return 0.5 * arma::sign(x) % (arma::abs(z) + z);
}

arma::vec soft_threshold(const arma::vec& x, const arma::vec& lambda) {
  arma::vec z = arma::abs(x) - lambda;
  return 0.5 * arma::sign(x) % (arma::abs(z) + z);
}


//' ADMM optimizer for generalized lasso problems
//' 
//' @param x (arma::mat) design matrix.
//' @param y (arma::vec) response vector.
//' @param w (arma::vec) weight vector.
//' @param nu (arma::vec) penalty weight vector.
//' @param pL2 (arma::vec) ridge penalty vector.
//' @param beta0 (arma::vec) initial values.
//' @param lambda (double) generalized lasso penalty parameter.
//' @param rho (double) augmented Lagrangian penalty parameter.
//' @param objtol (double) objective function tollerance parameter.
//' @param reltol (double) relative tollerance parameter.
//' @param abstol (double) absolute tollerance parameter.
//' @param maxiter (int) maximum numner of iterations.
//' 
//' @return A list containing the estimated coefficients and the optimization history.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List admm_genlasso_PQ(const arma::vec& y, 
                            const arma::mat& X, 
                            const arma::vec& w, 
                            const arma::vec& nu, 
                            const arma::vec& pL2, 
                            const arma::vec& beta0, 
                            const double& lambda, 
                            const double& rho, 
                            const bool& precondition,
                            const double& objtol,
                            const double& reltol, 
                            const double& abstol, 
                            const int& maxiter,
                            const bool& history) {
  
  // Data dimensions
  const int n = X.n_rows;
  const int p = X.n_cols;
  
  const double sqrtn = std::sqrt(static_cast<double>(n));
  const double sqrtp = std::sqrt(static_cast<double>(p));
  
  
  // Set the primal and dual variables and residuals
  arma::vec beta = beta0;
  arma::vec Xb = X * beta;
  arma::vec NXb = nu % Xb;
  arma::vec z = NXb;
  arma::vec z_old = z;
  arma::vec u(n, fill::zeros);
  arma::vec q(n, fill::zeros);
  arma::vec r(n, fill::zeros);
  arma::vec s(n, fill::zeros);
  arma::vec pw(n, fill::ones);
  arma::vec thr(n, fill::ones);
  
  // Set the preconditioning weights
  if(precondition){
    pw = 1 / (1e-6 + arma::abs(nu) % arma::sum(arma::abs(X), 1));
  }
  
  // Set the soft-thresholding operator parameters
  thr = (lambda/rho)/pw;
  
  // Set the primal and dual diagnostic variables
  double loss = 0.5*arma::accu(w % arma::square(y - Xb));
  double pridge = 0.5*arma::accu(pL2 % arma::square(beta));
  double plasso = lambda*arma::norm(NXb, 1); 
  
  double objval = loss + pridge + plasso;
  double objold = objval;
  double objmin = objval;
  double diffobj = .0;
  
  double norm_b = arma::norm(beta, 2);
  double norm_z = arma::norm(z, 2);
  double norm_u = .0;
  double norm_r = .0;
  double norm_s = .0;
  
  double eps_pri = .0;
  double eps_dual = .0;
  
  bool check_objval = false;
  bool check_norm_r = false;
  bool check_norm_s = false;
  
  // Compute the factorization
  arma::vec weights = w + rho * (nu % nu % pw);
  arma::mat cholQ;
  arma::mat XP;
  arma::mat XPXt;
  arma::mat XtWX;
  arma::vec wy = w % y;
  
  if(p>n){
    XP = X * diagmat(1/pL2);
    XPXt = X * XP.t();
    arma::chol(cholQ, XPXt + arma::diagmat(1/weights), "upper");
  }else{
    XtWX = X.t() * arma::diagmat(weights) * X;
    arma::chol(cholQ, XtWX + arma::diagmat(pL2), "upper");
  }
  
  // Set the optimization history
  arma::vec state(9);
  arma::mat trace(0, 9);
  
  if(history){
    state = arma::vec{.0, objval, norm_b, norm_z, norm_u, norm_r, norm_s, eps_pri, eps_dual};
    trace = arma::join_cols(trace, state.t());
  } 
  
  // ADMM optimization cycle
  int iter;
  for(iter=1; iter<maxiter; iter++){
    // Primal update: beta
    q = X.t() * (wy + rho * (pw % nu) % (z - u));
    if(p>n){
      beta = q/pL2 - XP.t() * chol_solve(cholQ, XP * q);
    }else{
      beta = chol_solve(cholQ, q);
    }
    
    // Primal update: eta
    Xb = X * beta;
    NXb = nu % Xb;
    
    // Primal update: z
    z_old = z;
    z = soft_threshold(NXb + u, thr);
    
    // Dual update: u
    u += - z + NXb;
    
    // Primal and dual residuals
    r = NXb - z;
    s = - (rho * pw) % (z - z_old);
    
    // Primal objective function
    loss = 0.5*arma::accu(w % arma::square(y - Xb));
    pridge = 0.5*arma::accu(pL2 % arma::square(beta));
    plasso = lambda*arma::norm(NXb, 1);
    
    objval = loss + pridge + plasso;
    objmin = std::min(objval, objold)/n;
    diffobj = std::abs(objval - objold)/n;
    
    // Convergence indicators
    norm_b = arma::norm(beta, 2);
    norm_z = arma::norm(z, 2);
    norm_u = arma::norm(u, 2);
    norm_r = arma::norm(r, 2);
    norm_s = arma::norm(s, 2);
    eps_pri = sqrtp*abstol + reltol*std::max(norm_b, norm_z);
    eps_dual = sqrtn*abstol + reltol*rho*norm_u;
    
    // Store the optimization state
    if(history){
      state = arma::vec{double(iter), objval, norm_b, norm_z, norm_u, norm_r, norm_s, eps_pri, eps_dual};
      trace = arma::join_cols(trace, state.t());
    }
    
    // Convergence checks
    check_objval = (diffobj < n*objtol) & (objval <= objmin);
    check_norm_r = norm_r < eps_pri;
    check_norm_s = norm_s < eps_dual;
    if(check_objval & check_norm_r & check_norm_s){
      break;
    }
  }
  
  // Output list
  Rcpp::List output;
  output["beta"] = beta;
  output["z"] = z;
  output["u"] = u;
  output["niter"] = iter;
  if(history){
    output["trace"] = trace;
  }
  
  return output;
}
  



