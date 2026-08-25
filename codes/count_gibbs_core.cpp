// mncount_gibbs_core.cpp
//
// Count-model extensions for the Pólya-Gamma Gibbs sampler.
// Compile alongside (or append to) mnlogit_gibbs_core.cpp.
//
// New exports vs. MNL core:
//   compute_utilities_count_cpp  — linear predictor only (no softmax / c_j)
//   gibbs_step_count_cpp         — beta update without competitor term
//   compute_loglik_nb_cpp        — Negative-Binomial log-likelihood
//   compute_loglik_poisson_cpp   — Poisson log-likelihood (for diagnostics)
//
// All other functions (chol_sample_precision_cpp, weighted_crossprod,
// update_re_precision, sample_tmvn_precision_gibbs_cpp, …) are shared
// with the MNL core and are NOT duplicated here.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Forward-declare shared helpers (defined in mnlogit_gibbs_core.cpp)
arma::vec chol_sample_precision_cpp(const arma::mat &P, const arma::vec &Pb);

// ==========================================================================
// 1. Linear predictor only (count model: no softmax, no c_j)
//
//    U[i, ip] = X[i,] %*% beta[, ip]  +  f_bart[i, ip]  +  offset[i, ip]
//
//    All p columns are returned directly (no baseline removal).
// ==========================================================================
// [[Rcpp::export]]
arma::mat compute_utilities_count_cpp(
    const arma::mat &X,
    const arma::mat &beta,
    const arma::mat &f_bart,
    const arma::mat &offset)
{
  int n = X.n_rows;
  int p = beta.n_cols;
  arma::mat U(n, p);
  for (int ip = 0; ip < p; ip++) {
    U.col(ip) = X * beta.col(ip) + f_bart.col(ip) + offset.col(ip);
    // Wider clamp than MNL (counts can legitimately need larger |eta|)
    U.col(ip) = arma::clamp(U.col(ip), -30.0, 30.0);
  }
  return U;
}

// ==========================================================================
// 2. Full inner Gibbs step for POOLED count model
//
//    KEY DIFFERENCE FROM MNL (gibbs_step_pooled):
//      • No c_j term.  target = kappa[, ip] (– omega * f_bart if use_bart).
//      • All p columns are active; no pp index mapping needed.
//
//    Arguments:
//      kappa_mat  n x p  matrix  =  (Y - r_j) / 2   (updated per iteration)
//      omega      n x p  Pólya-Gamma augmentation weights
//      f_bart     n x p  current BART function values (zeros if !use_bart)
// ==========================================================================
// [[Rcpp::export]]
arma::mat gibbs_step_count_cpp(
    const arma::mat &X,
    const arma::mat &Xt,
    const arma::mat &kappa_mat,
    const arma::mat &omega,
    const arma::mat &prior_P,
    const arma::mat &hs_prec,
    const arma::mat &prior_Pb,
    const arma::mat &f_bart,
    bool use_bart)
{
  int k = Xt.n_rows;
  int p = omega.n_cols;
  arma::mat beta_new(k, p);

  for (int ip = 0; ip < p; ip++) {
    arma::vec om_p = omega.col(ip);

    // Precision: prior  +  X' diag(omega) X  +  Horseshoe correction
    arma::mat P = prior_P + Xt * (X.each_col() % om_p);
    P.diag()   += hs_prec.col(ip);

    // Right-hand side: kappa  (minus omega * f_bart if BART component present)
    // Derivation: working response is (kappa/omega), precision-weighted RHS is
    //   X' diag(omega) * (kappa/omega) = X' kappa
    // With BART: working response is (kappa - omega*f_bart)/omega,
    //   so RHS = X' (kappa - omega*f_bart)
    arma::vec target = kappa_mat.col(ip);
    if (use_bart)
      target -= om_p % f_bart.col(ip);

    arma::vec Pb = prior_Pb.col(ip) + Xt * target;
    beta_new.col(ip) = chol_sample_precision_cpp(P, Pb);
  }
  return beta_new;
}

// ==========================================================================
// 3. Negative-Binomial log-likelihood (logit parameterisation)
//
//    Model:  Y_{ij} ~ NB(r_j, sigmoid(eta_{ij}))
//    Parameterisation:
//      log p(y|r,eta) = lgamma(y+r) − lgamma(r) − lgamma(y+1)
//                       + y*eta − (y+r)*log1p(exp(eta))
//
//    r_vec: length-p vector of dispersion parameters.
//
//    Returns:
//      total     : scalar total log-likelihood (sum over i and j)
//      pointwise : n-vector of row-summed log-likelihoods (for LOO-CV)
// ==========================================================================
// [[Rcpp::export]]
List compute_loglik_nb_cpp(
    const arma::mat &eta,   // n x p  linear predictor
    const arma::mat &Y,     // n x p  count matrix
    const arma::vec &r_vec) // length p
{
  int n = eta.n_rows;
  int p = eta.n_cols;
  arma::vec ll_pw(n, arma::fill::zeros);

  for (int ip = 0; ip < p; ip++) {
    double r = r_vec(ip);
    for (int i = 0; i < n; i++) {
      double y   = Y(i, ip);
      double e   = eta(i, ip);
      // log p(y | r, eta)
      double lp  = R::lgammafn(y + r) - R::lgammafn(r) - R::lgammafn(y + 1.0)
                   + y * e - (y + r) * std::log1p(std::exp(e));
      if (std::isfinite(lp)) ll_pw(i) += lp;
    }
  }
  return List::create(Named("total")     = arma::accu(ll_pw),
                      Named("pointwise") = ll_pw);
}

// ==========================================================================
// 4. Poisson log-likelihood (log-link)
//
//    Model:  Y_{ij} ~ Poisson(exp(eta_{ij}))
//    log p(y|eta) = y*eta − exp(eta) − lgamma(y+1)
//
//    Useful for diagnostics when family = "poisson" (r large → NB → Poisson).
// ==========================================================================
// [[Rcpp::export]]
List compute_loglik_poisson_cpp(
    const arma::mat &eta,
    const arma::mat &Y)
{
  int n = eta.n_rows;
  int p = eta.n_cols;
  arma::vec ll_pw(n, arma::fill::zeros);

  for (int ip = 0; ip < p; ip++) {
    for (int i = 0; i < n; i++) {
      double y  = Y(i, ip);
      double e  = eta(i, ip);
      double lp = y * e - std::exp(e) - R::lgammafn(y + 1.0);
      if (std::isfinite(lp)) ll_pw(i) += lp;
    }
  }
  return List::create(Named("total")     = arma::accu(ll_pw),
                      Named("pointwise") = ll_pw);
}

// ==========================================================================
// 5. Predicted mean (count scale) from posterior draws
//
//    For NB with logit parameterisation:
//      mu_{ij} = r_j * exp(eta_{ij})    (mean of the NB distribution)
//      E[Y] = r * p / (1-p) = r * sigmoid(eta) / (1 - sigmoid(eta))
//           = r * exp(eta)
//
//    eta_draws: n x p x S array (S posterior samples), passed as n x (p*S)
//               column-major; caller reshapes in R.
//    r_draws  : p x S matrix of dispersion draws.
//
//    Returns n x p matrix of posterior mean predictions.
// ==========================================================================
// [[Rcpp::export]]
arma::mat posterior_mean_nb_cpp(
    const arma::mat &eta,   // n x p  (single draw or posterior mean eta)
    const arma::vec &r_vec) // length p
{
  int n = eta.n_rows;
  int p = eta.n_cols;
  arma::mat mu(n, p);
  for (int ip = 0; ip < p; ip++) {
    double r = r_vec(ip);
    for (int i = 0; i < n; i++) {
      mu(i, ip) = r * std::exp(eta(i, ip));
    }
  }
  return mu;
}
