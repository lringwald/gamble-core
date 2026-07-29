// mnlogit_gibbs_core.cpp
// RcppArmadillo implementation of the inner Gibbs loop for the MNL sampler.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ==========================================================================
// 1. Cholesky sample from precision parameterization
// ==========================================================================
// [[Rcpp::export]]
arma::vec chol_sample_precision_cpp(const arma::mat &P, const arma::vec &Pb) {
  int k = P.n_rows;
  arma::mat R;
  bool ok = arma::chol(R, P);
  if (!ok) {
    arma::mat P_jit = P + 1e-6 * arma::eye(k, k);
    ok = arma::chol(R, P_jit);
    if (!ok) {
      arma::mat V = arma::pinv(P);
      arma::vec mu = V * Pb;
      return mu + arma::chol(arma::symmatu(V), "lower") * arma::randn(k);
    }
  }
  arma::vec mu =
      arma::solve(arma::trimatu(R), arma::solve(arma::trimatl(R.t()), Pb));
  arma::vec z = arma::randn(k);
  return mu + arma::solve(arma::trimatu(R), z);
}

// ==========================================================================
// 2. Cholesky: return mean and covariance
// ==========================================================================
// [[Rcpp::export]]
List chol_mean_and_var_cpp(const arma::mat &P, const arma::vec &Pb) {
  int k = P.n_rows;
  arma::mat R;
  bool ok = arma::chol(R, P);
  if (!ok) {
    arma::mat P_jit = P + 1e-6 * arma::eye(k, k);
    ok = arma::chol(R, P_jit);
    if (!ok) {
      arma::mat V = arma::pinv(P);
      V = 0.5 * (V + V.t());
      arma::vec mu = V * Pb;
      return List::create(Named("mu") = mu, Named("V") = V);
    }
  }
  arma::vec mu =
      arma::solve(arma::trimatu(R), arma::solve(arma::trimatl(R.t()), Pb));
  arma::mat R_inv = arma::solve(arma::trimatu(R), arma::eye(k, k));
  arma::mat V = R_inv * R_inv.t();
  return List::create(Named("mu") = mu, Named("V") = V);
}

// ==========================================================================
// 3. Build precision matrix
// ==========================================================================
// [[Rcpp::export]]
arma::mat weighted_crossprod(const arma::mat &Xt, const arma::mat &X,
                             const arma::vec &w) {
  return Xt * (X.each_col() % w);
}

// ==========================================================================
// 4. Utility update + log-sum-exp + c_j computation
// ==========================================================================
// [[Rcpp::export]]
List update_utilities_and_cj(const arma::mat &X, const arma::mat &beta,
                             const arma::mat &f_bart, const arma::ivec &pp,
                             int baseline, int p_all, bool use_bart) {
  int n = X.n_rows;
  int p = beta.n_cols;
  arma::mat U(n, p_all, arma::fill::zeros);
  for (int ip = 0; ip < p; ip++) {
    int j = pp(ip) - 1;
    U.col(j) = X * beta.col(ip);
    if (use_bart)
      U.col(j) += f_bart.col(ip);
    U.col(j) = arma::clamp(U.col(j), -20.0, 20.0);
  }
  arma::vec m_U = arma::max(U, 1);
  arma::mat U_shifted = U.each_col() - m_U;
  arma::mat exp_U = arma::exp(U_shifted);
  arma::vec sum_exp = arma::sum(exp_U, 1);
  arma::mat c_j_mat(n, p);
  for (int ip = 0; ip < p; ip++) {
    int j = pp(ip) - 1;
    arma::vec remainder = sum_exp - exp_U.col(j);
    remainder = arma::clamp(remainder, 1e-10, datum::inf);
    c_j_mat.col(ip) = m_U + arma::log(remainder);
  }
  return List::create(Named("U") = U, Named("c_j_mat") = c_j_mat);
}

// [[Rcpp::export]]
List update_utilities_and_cj_re(const List &Xm_list, const List &idx_list,
                                int n, const arma::cube &beta_c,
                                const arma::mat &f_bart, const arma::ivec &pp,
                                int baseline, int p_all, bool use_bart) {
  int p = beta_c.n_cols;
  int n_groups = beta_c.n_slices;
  arma::mat U(n, p_all, arma::fill::zeros);

  // Per-group linear predictor via one BLAS-3 GEMM on the design matrix the
  // caller already pre-gathered once (Xm_list[m]), replacing the scalar
  // for-i/for-ip/for-v cube loop AND avoiding any per-call row gather.
  for (int m = 0; m < n_groups; m++) {
    arma::uvec idx = as<arma::uvec>(idx_list[m]) - 1; // 1-based -> 0-based
    if (idx.n_elem == 0)
      continue;
    arma::mat Xm = as<arma::mat>(Xm_list[m]);
    arma::mat Eta = Xm * beta_c.slice(m); // n_m x p
    if (use_bart)
      Eta += f_bart.rows(idx);
    Eta = arma::clamp(Eta, -20.0, 20.0);
    for (int ip = 0; ip < p; ip++) {
      int j = pp(ip) - 1;
      U(idx, arma::uvec({(unsigned int)j})) = Eta.col(ip);
    }
  }

  arma::vec m_U = arma::max(U, 1);
  arma::mat U_shifted = U.each_col() - m_U;
  arma::mat exp_U = arma::exp(U_shifted);
  arma::vec sum_exp = arma::sum(exp_U, 1);
  arma::mat c_j_mat(n, p);
  for (int ip = 0; ip < p; ip++) {
    int j = pp(ip) - 1;
    arma::vec remainder = sum_exp - exp_U.col(j);
    remainder = arma::clamp(remainder, 1e-10, datum::inf);
    c_j_mat.col(ip) = m_U + arma::log(remainder);
  }

  return List::create(Named("U") = U, Named("c_j_mat") = c_j_mat);
}

// ==========================================================================
// 5. Full inner Gibbs step for POOLED model
// ==========================================================================
// [[Rcpp::export]]
arma::mat gibbs_step_pooled(const arma::mat &X, const arma::mat &Xt,
                            const arma::mat &kappa_w, const arma::mat &omega,
                            const arma::mat &c_j_mat, const arma::mat &prior_P,
                            const arma::mat &hs_prec, const arma::mat &prior_Pb,
                            const arma::ivec &pp, const arma::mat &f_bart,
                            bool use_bart) {
  int k = Xt.n_rows;
  int p = omega.n_cols;
  arma::mat beta_new(k, p);
  for (int ip = 0; ip < p; ip++) {
    int j = pp(ip) - 1;
    arma::vec om_p = omega.col(ip);
    arma::vec c_j = c_j_mat.col(ip);
    arma::mat P = prior_P + Xt * (X.each_col() % om_p);
    P.diag() += hs_prec.col(ip);
    arma::vec target = kappa_w.col(j) + om_p % c_j;
    if (use_bart)
      target -= om_p % f_bart.col(ip);
    arma::vec Pb = prior_Pb.col(ip) + Xt * target;
    beta_new.col(ip) = chol_sample_precision_cpp(P, Pb);
  }
  return beta_new;
}

// ==========================================================================
// 6. Full inner Gibbs step for RE model (Mixed Effects Support)
// ==========================================================================
// [[Rcpp::export]]
List gibbs_step_re(const arma::mat &X, const arma::mat &Xt,
                   const arma::mat &kappa_w, const arma::mat &omega,
                   const arma::mat &c_j_mat, const arma::mat &prior_P,
                   const arma::mat &hs_prec_mu, const arma::mat &prior_Pb,
                   const arma::mat &mu_pooled, const arma::mat &prec_pooled,
                   const arma::ivec &pp, const List &idx_list,
                   const List &Xm_list, const List &Xmt_list, int n_groups,
                   const arma::mat &f_bart, bool use_bart,
                   const arma::uvec &re_idx_0, const arma::mat &re_mask,
                   const arma::mat &y_mask, const arma::uvec &is_intercept) {
  int k = X.n_cols;
  int p = omega.n_cols;
  arma::cube beta_c(k, p, n_groups, arma::fill::zeros);
  arma::mat mu_new(k, p);

  // Convert re_idx to a logical mask for fast lookup
  arma::uvec is_re(k, arma::fill::zeros);
  for (unsigned int i = 0; i < re_idx_0.n_elem; i++)
    is_re(re_idx_0(i)) = 1;

  for (int ip = 0; ip < p; ip++) {
    int j = pp(ip) - 1;
    arma::vec om_p = omega.col(ip);
    arma::vec c_j = c_j_mat.col(ip);
    arma::vec diag_p = prec_pooled.col(ip);

    // Accumulators for residual-based pooled update
    arma::mat P_fixed_sum(k, k, arma::fill::zeros);
    arma::vec resid_Pb_sum(k, arma::fill::zeros);

    for (int m = 0; m < n_groups; m++) {
      arma::uvec idx = as<arma::uvec>(idx_list[m]) - 1;
      if (idx.n_elem == 0)
        continue;

      arma::mat Xm = as<arma::mat>(Xm_list[m]);
      arma::mat Xmt = as<arma::mat>(Xmt_list[m]);
      arma::vec om_m = om_p(idx);

      arma::vec target_m =
          kappa_w(idx, arma::uvec({(unsigned int)j})) + om_m % c_j(idx);
      if (use_bart)
        target_m -= om_m % f_bart(idx, arma::uvec({(unsigned int)ip}));

      // A. Sample hierarchical coefficients beta_g
      // Note: For FIXED indices, diag_p is very large (1e12), so beta_g will be
      // pinned to mu_pooled. This is okay for the group step.
      arma::mat P_m = arma::diagmat(diag_p) + Xmt * (Xm.each_col() % om_m);
      arma::vec Pb_m = diag_p % mu_pooled.col(ip) + Xmt * target_m;
      beta_c.slice(m).col(ip) = chol_sample_precision_cpp(P_m, Pb_m);

      for (int v = 0; v < k; v++) {
        if (!is_re(v) || re_mask(v, m) == 0.0 || (is_intercept(v) == 0 && y_mask(ip, m) == 0.0))
          beta_c(v, ip, m) = mu_pooled(v, ip);
      }
      // B. Accumulate for the global mu update
      P_fixed_sum += Xmt * (Xm.each_col() % om_m);

      // Residual-based RHS: target - X_re * delta_g, where delta_g = beta_g -
      // mu This allows the global mu to capture the pooled signal correctly.
      arma::vec delta_g = beta_c.slice(m).col(ip) - mu_pooled.col(ip);
      arma::vec resid_target = target_m - (Xm * delta_g) % om_m;
      resid_Pb_sum += Xmt * resid_target;


    }

    // --- Sample mu_pooled (Global parameters: Fixed Effects + Hierarchical
    // Means) ---
    arma::mat P_mu = prior_P + P_fixed_sum;
    arma::vec Pb_mu = prior_Pb.col(ip) + resid_Pb_sum;
    P_mu.diag() += hs_prec_mu.col(ip);

    mu_new.col(ip) = chol_sample_precision_cpp(P_mu, Pb_mu);

    // --- Final Pass: Ensure fixed beta_g are pinned to mu_new ---
    for (int m = 0; m < n_groups; m++) {
      for (int v = 0; v < k; v++) {
        if (!is_re(v) || re_mask(v, m) == 0.0 || (is_intercept(v) == 0 && y_mask(ip, m) == 0.0))
          beta_c(v, ip, m) = mu_new(v, ip);
      }
    }
  }
  return List::create(Named("beta_c") = beta_c, Named("mu") = mu_new);
}

// ==========================================================================
// 7. Sigma (group-level precision) update
// ==========================================================================
// [[Rcpp::export]]
List update_re_precision(const arma::cube &beta_c, const arma::mat &mu_pooled,
                         const arma::uvec &re_idx, int n_groups,
                         const arma::cube &re_mask, const arma::mat &y_mask,
                         const arma::uvec &is_intercept, double prior_a_re = 0.5,
                         double prior_b_re = 0.5) {
  int k = beta_c.n_rows;
  int p = beta_c.n_cols;
  arma::mat prec_out(k, p);
  prec_out.fill(1e12);
  arma::mat sigma_out(k, p);
  sigma_out.fill(1e-6);
  for (int ip = 0; ip < p; ip++) {
    for (unsigned int ri = 0; ri < re_idx.n_elem; ri++) {
      int v = re_idx(ri);
      double ss = 0.0;
      double n_active = 0.0;
      for (int m = 0; m < n_groups; m++) {
        if (re_mask(v, ip, m) > 0.5 && (is_intercept(v) == 1 || y_mask(ip, m) > 0.5)) {
          double diff = beta_c(v, ip, m) - mu_pooled(v, ip);
          ss += diff * diff;
          n_active += 1.0;
        }
      }
      if (n_active < 1.0) {
        prec_out(v, ip) = 1.0;
        sigma_out(v, ip) = 1.0;
        continue;
      }
      double rate = prior_b_re + ss / 2.0;
      if (!std::isfinite(rate) || rate <= 0)
        rate = prior_b_re;
      double shape = prior_a_re + n_active / 2.0;
      double prec = R::rgamma(shape, 1.0 / rate);
      prec = std::min(std::max(prec, 1e-4), 1e6);
      prec_out(v, ip) = prec;
      sigma_out(v, ip) = 1.0 / std::sqrt(prec + 1e-12);
    }
  }
  return List::create(Named("prec") = prec_out, Named("sigma") = sigma_out);
}

// ==========================================================================
// 8. Log-likelihood computation
// ==========================================================================
// [[Rcpp::export]]
List compute_loglik(const arma::mat &U, const arma::mat &Y,
                    const arma::vec &y_weight_vec) {
  arma::vec m_U = arma::max(U, 1);
  arma::mat U_shifted = U.each_col() - m_U;
  arma::vec lse = m_U + arma::log(arma::sum(arma::exp(U_shifted), 1));
  arma::mat log_probs = U.each_col() - lse;
  arma::vec ll_i = arma::sum(Y % log_probs, 1);

  arma::vec ll_weighted;
  if (y_weight_vec.n_elem == 1) {
    ll_weighted = ll_i * y_weight_vec[0];
  } else {
    ll_weighted = ll_i % y_weight_vec;
  }

  double total_ll = arma::accu(ll_weighted);
  return List::create(Named("total") = total_ll,
                      Named("pointwise") = ll_weighted);
}

// ==========================================================================
// 9. BART Prediction Engine (Hardened)
// ==========================================================================
// [[Rcpp::export]]
arma::vec predict_slim_bart_cpp(const arma::mat &X, const DataFrame &tree_df) {
  int n_obs = X.n_rows;
  IntegerVector tree_ids = tree_df["tree"];
  IntegerVector vars = tree_df["var"];
  NumericVector vals = tree_df["value"];
  int total_rows = tree_df.nrows();

  std::vector<int> tree_starts;
  if (total_rows > 0) {
    tree_starts.push_back(0);
    for (int i = 1; i < total_rows; i++) {
      if (tree_ids[i] != tree_ids[i - 1])
        tree_starts.push_back(i);
    }
  }

  std::vector<int> right_children(total_rows, -1);
  auto get_subtree_size = [&](int start_idx, auto &self) -> int {
    if (start_idx >= total_rows)
      return 0;
    if (vars[start_idx] == -1)
      return 1;
    int left_child = start_idx + 1;
    int left_size = self(left_child, self);
    int right_child = left_child + left_size;
    right_children[start_idx] = right_child;
    int right_size = self(right_child, self);
    return 1 + left_size + right_size;
  };

  for (int start_node : tree_starts)
    get_subtree_size(start_node, get_subtree_size);

  arma::vec out(n_obs, arma::fill::zeros);
  for (int i = 0; i < n_obs; i++) {
    for (int start_node : tree_starts) {
      int curr = start_node;
      while (vars[curr] != -1) {
        int v = vars[curr] - 1;
        if (X(i, v) <= vals[curr])
          curr++;
        else
          curr = right_children[curr];
        if (curr == -1 || curr >= total_rows)
          break;
      }
      if (curr != -1 && curr < total_rows)
        out[i] += vals[curr];
    }
  }
  return out;
}
// ==========================================================================
// 10. Univariate truncated normal via inverse-CDF
//     Safe against numerical edge cases at the tails
// ==========================================================================
inline double rtruncnorm_scalar(double mean, double sd, double lo, double hi) {
  // Standardize the bounds
  double a = (lo - mean) / sd;
  double b = (hi - mean) / sd;

  // Helper: Robert (1995) optimal exponential rejection for the tail
  auto rnorm_tail = [](double a) {
    double alpha = (a + std::sqrt(a * a + 4.0)) / 2.0;
    while (true) {
      double z = a - std::log(1.0 - R::unif_rand()) / alpha;
      double diff = z - alpha;
      double rho = std::exp(-0.5 * diff * diff);
      if (R::unif_rand() < rho)
        return z;
    }
  };

  double res;

  if (a > 0.5) {
    // Case 1: Right Tail
    if (b == R_PosInf) {
      res = rnorm_tail(a);
    } else {
      // Finite interval in tail: Rejection from Uniform
      while (true) {
        res = a + R::unif_rand() * (b - a);
        if (R::unif_rand() < std::exp(0.5 * (a * a - res * res)))
          break;
      }
    }
  } else if (b < -0.5) {
    // Case 2: Left Tail (Symmetry)
    if (a == R_NegInf) {
      res = -rnorm_tail(-b);
    } else {
      while (true) {
        res = b - R::unif_rand() * (b - a);
        if (R::unif_rand() < std::exp(0.5 * (b * b - res * res)))
          break;
      }
    }
  } else {
    // Case 3: Central Region
    // Optimization: If interval is very narrow, use Uniform rejection
    // instead of Normal rejection to avoid infinite loops.
    if (std::abs(b - a) < 0.2) {
      while (true) {
        res = a + R::unif_rand() * (b - a);
        if (R::unif_rand() < std::exp(-0.5 * res * res))
          break;
      }
    } else {
      while (true) {
        res = R::rnorm(0, 1);
        if (res >= a && res <= b)
          break;
      }
    }
  }

  return mean + sd * res;
}

// ==========================================================================
// 11. Coordinate Gibbs sampler for truncated MVN (precision parameterization)
//     P   : k x k precision matrix
//     Pb  : k-vector = P %*% mu  (already computed upstream — pass it in)
//     lo  : k-vector lower bounds (-Inf for unconstrained)
//     hi  : k-vector upper bounds (+Inf for unconstrained)
//     init: k-vector starting value (previous draw)
//     n_steps: number of full coordinate sweeps (1-3 usually sufficient)
// ==========================================================================
// [[Rcpp::export]]
arma::vec
sample_tmvn_precision_gibbs_cpp(const arma::mat &P, const arma::vec &Pb,
                                const arma::vec &lo, const arma::vec &hi,
                                const arma::vec &init, int n_steps = 2) {
  int k = P.n_rows;

  // Compute mean once: mu = P^{-1} Pb via Cholesky
  arma::mat R;
  arma::vec mu;
  bool ok = arma::chol(R, P);
  if (ok) {
    mu = arma::solve(arma::trimatu(R), arma::solve(arma::trimatl(R.t()), Pb));
  } else {
    // Fallback: pinv (should rarely trigger)
    mu = arma::pinv(P) * Pb;
  }

  arma::vec x = init;

  // Ensure init is feasible
  for (int j = 0; j < k; j++) {
    x(j) = std::max(lo(j), std::min(hi(j), x(j)));
  }

  GetRNGstate();

  for (int step = 0; step < n_steps; step++) {
    for (int j = 0; j < k; j++) {
      double Pjj = P(j, j);

      // Conditional mean: mu_j - (1/Pjj) * P[j,-j] . (x[-j] - mu[-j])
      // Computed as a dot product, skipping index j
      double correction = 0.0;
      for (int l = 0; l < k; l++) {
        if (l != j)
          correction += P(j, l) * (x(l) - mu(l));
      }
      double cond_mean = mu(j) - correction / Pjj;
      double cond_sd = 1.0 / std::sqrt(Pjj);

      x(j) = rtruncnorm_scalar(cond_mean, cond_sd, lo(j), hi(j));
    }
  }

  PutRNGstate();
  return x;
}
// ==========================================================================
// 12. Inner Gibbs step for RE model — Non-Centered Parameterization
//
//  Model:  beta_c[, m] = mu + sigma * z_c[, m],  z_c[, m] ~ N(0, I)
//
//  Sampling order per category ip:
//    A. For each group m: draw z_c[, m] from its conditional
//    B. Draw mu from its conditional (z-residualized sufficient stats)
//    C. Recover beta_c = mu + sigma * z_c and pin fixed effects
//
//  Arguments:
// ==========================================================================
// [[Rcpp::export]]
List gibbs_step_re_ncp(SEXP X_s, SEXP Xt_s, SEXP kappa_w_s, SEXP omega_s,
                       SEXP c_j_mat_s, SEXP prior_P_s, SEXP hs_prec_mu_s,
                       SEXP prior_Pb_s, SEXP mu_pooled_s, SEXP sigma_mat_s,
                       SEXP z_c_in_s, SEXP pp_s, SEXP idx_list_s,
                       SEXP Xm_list_s, SEXP Xmt_list_s, int n_groups,
                       SEXP f_bart_s, bool use_bart, SEXP re_idx_0_s,
                       SEXP re_mask_s, SEXP y_mask_s, SEXP is_intercept_s,
                       SEXP re_support_s, SEXP a_re_s, bool re_asis) {
  NumericMatrix X_m = as<NumericMatrix>(X_s);
  NumericMatrix Xt_m = as<NumericMatrix>(Xt_s);
  NumericMatrix kappa_w_m = as<NumericMatrix>(kappa_w_s);
  NumericMatrix omega_m = as<NumericMatrix>(omega_s);
  NumericMatrix c_j_mat_m = as<NumericMatrix>(c_j_mat_s);
  NumericMatrix prior_P_m = as<NumericMatrix>(prior_P_s);
  NumericMatrix hs_prec_mu_m = as<NumericMatrix>(hs_prec_mu_s);
  NumericMatrix prior_Pb_m = as<NumericMatrix>(prior_Pb_s);
  NumericMatrix mu_pooled_m = as<NumericMatrix>(mu_pooled_s);
  NumericMatrix sigma_mat_m = as<NumericMatrix>(sigma_mat_s);
  NumericVector z_c_in_v = as<NumericVector>(z_c_in_s);
  IntegerVector pp_v = as<IntegerVector>(pp_s);
  List idx_list = as<List>(idx_list_s);
  List Xm_list = as<List>(Xm_list_s);
  List Xmt_list = as<List>(Xmt_list_s);
  NumericMatrix f_bart_m = as<NumericMatrix>(f_bart_s);
  IntegerVector re_idx_0_v = as<IntegerVector>(re_idx_0_s);
  NumericMatrix re_mask_m = as<NumericMatrix>(re_mask_s);
  NumericMatrix y_mask_m = as<NumericMatrix>(y_mask_s);
  LogicalVector is_intercept_v = as<LogicalVector>(is_intercept_s);

  // 1. Bridge to Armadillo
  arma::mat X(X_m.begin(), X_m.nrow(), X_m.ncol(), false);
  arma::mat Xt(Xt_m.begin(), Xt_m.nrow(), Xt_m.ncol(), false);
  arma::mat kappa_w(kappa_w_m.begin(), kappa_w_m.nrow(), kappa_w_m.ncol(),
                    false);
  arma::mat omega(omega_m.begin(), omega_m.nrow(), omega_m.ncol(), false);
  arma::mat c_j_mat(c_j_mat_m.begin(), c_j_mat_m.nrow(), c_j_mat_m.ncol(),
                    false);
  arma::mat prior_P(prior_P_m.begin(), prior_P_m.nrow(), prior_P_m.ncol(),
                    false);
  arma::mat hs_prec_mu(hs_prec_mu_m.begin(), hs_prec_mu_m.nrow(),
                       hs_prec_mu_m.ncol(), false);
  arma::mat prior_Pb(prior_Pb_m.begin(), prior_Pb_m.nrow(), prior_Pb_m.ncol(),
                      false);
  arma::mat mu_pooled(mu_pooled_m.begin(), mu_pooled_m.nrow(),
                      mu_pooled_m.ncol(), false);
  arma::mat sigma_mat(sigma_mat_m.begin(), sigma_mat_m.nrow(),
                      sigma_mat_m.ncol(), false);
  arma::mat re_mask(re_mask_m.begin(), re_mask_m.nrow(), re_mask_m.ncol(),
                    false);
  arma::mat y_mask(y_mask_m.begin(), y_mask_m.nrow(), y_mask_m.ncol(), false);
  arma::uvec is_intercept = as<arma::uvec>(is_intercept_v);
  // Support-aware RE prior: per (covariate x group) SD multiplier in (0,1] from the
  // within-group participation ratio (shrinks info-sparse RE slopes). All-ones = off.
  NumericMatrix re_support_m = as<NumericMatrix>(re_support_s);
  arma::mat re_support(re_support_m.begin(), re_support_m.nrow(),
                       re_support_m.ncol(), false);
  // Half-Cauchy auxiliary a (sigma^2 ~ IG(1/2, 1/a)) for the RE-ASIS scale interweave.
  NumericMatrix a_re_m = as<NumericMatrix>(a_re_s);
  arma::mat a_re(a_re_m.begin(), a_re_m.nrow(), a_re_m.ncol(), false);

  IntegerVector z_dims = z_c_in_v.attr("dim");
  arma::cube z_c_in(z_c_in_v.begin(), z_dims[0], z_dims[1], z_dims[2], false);

  arma::ivec pp(pp_v.begin(), pp_v.size(), false);

  arma::mat f_bart(f_bart_m.begin(), f_bart_m.nrow(), f_bart_m.ncol(), false);
  arma::uvec re_idx_0 = as<arma::uvec>(re_idx_0_v);

  int k = X.n_cols;
  int p = omega.n_cols;
  int n = X.n_rows;

  // 2. Prepare output storage
  arma::cube z_c_out(k, p, n_groups, arma::fill::zeros);
  arma::cube beta_c(k, p, n_groups, arma::fill::zeros);
  arma::mat mu_new(k, p);

  // 3. Pre-cast lists to std::vectors for stability
  std::vector<IntegerVector> idx_vec(n_groups);
  std::vector<NumericMatrix> Xm_vec(n_groups);
  std::vector<NumericMatrix> Xmt_vec(n_groups);
  for (int m = 0; m < n_groups; m++) {
    idx_vec[m] = as<IntegerVector>(idx_list[m]);
    Xm_vec[m] = as<NumericMatrix>(Xm_list[m]);
    Xmt_vec[m] = as<NumericMatrix>(Xmt_list[m]);
  }

  // 4. Boolean mask for random effects
  arma::uvec is_re(k, arma::fill::zeros);
  for (unsigned int i = 0; i < re_idx_0.n_elem; i++)
    is_re(re_idx_0(i)) = 1;

  // 5. Gibbs Sampling
  for (int ip = 0; ip < p; ip++) {
    int j = pp(ip) - 1;

    arma::vec om_p = omega.col(ip);
    arma::vec c_j = c_j_mat.col(ip);
    arma::vec sig = sigma_mat.col(ip);
    arma::vec mu_ip = mu_pooled.col(ip);

    arma::mat P_mu_acc(k, k, arma::fill::zeros);
    arma::vec Pb_mu_acc(k, arma::fill::zeros);

    // Pre-declare variables outside the group loop for memory buffer reuse
    arma::vec sig_m;
    arma::mat G_m;
    arma::vec sig_red;
    arma::vec bdev;
    arma::vec r_z;
    arma::vec z_m;
    arma::uvec active_re;
    arma::mat P_z_red;
    arma::vec Pb_z_red;
    arma::vec z_red;
    arma::vec r_mu;

    for (int m = 0; m < n_groups; m++) {
      IntegerVector iv = idx_vec[m];
      arma::uvec idx = as<arma::uvec>(iv) - 1;
      if (idx.n_elem == 0) {
        z_c_out.slice(m).col(ip) = z_c_in.slice(m).col(ip);
        continue;
      }

      // Manual bridging from pre-casted NumericMatrix (no memory allocation)
      NumericMatrix Xm_m = Xm_vec[m];
      NumericMatrix Xmt_m = Xmt_vec[m];
      arma::mat Xm(Xm_m.begin(), Xm_m.nrow(), Xm_m.ncol(), false);
      arma::mat Xmt(Xmt_m.begin(), Xmt_m.nrow(), Xmt_m.ncol(), false);

      arma::vec om_m = om_p(idx);
      
      // Localize scales with structural mask + support-aware shrinkage
      sig_m = sig % re_mask.col(m) % re_support.col(m);
      if (y_mask(ip, m) < 0.5) {
        for (int v = 0; v < k; v++) {
          if (is_intercept(v) == 0) {
            sig_m(v) = 0.0;
          }
        }
      }

      // Single weighted Gram per group: G_m = Xm' diag(om) Xm.
      // Reused for BOTH the z-precision and the mu accumulation, replacing the
      // two separate crossprods in the original code.
      G_m = Xmt * (Xm.each_col() % om_m);

      r_z = kappa_w(idx, arma::uvec({(unsigned int)j})) +
            om_m % c_j(idx) - om_m % (Xm * mu_ip);
      if (use_bart)
        r_z -= om_m % f_bart(idx, arma::uvec({(unsigned int)ip}));

      z_m = arma::zeros<arma::vec>(k);

      // Identify active RE indices for group m
      std::vector<int> active_idx;
      for (int v = 0; v < k; v++) {
        bool active = is_re(v) && re_mask(v, m) > 0.5;
        if (active && (is_intercept(v) == 1 || y_mask(ip, m) > 0.5)) {
          active_idx.push_back(v);
        }
      }

      if (active_idx.size() > 0) {
        active_re.set_size(active_idx.size());
        for (size_t i = 0; i < active_idx.size(); i++) {
          active_re(i) = active_idx[i];
        }

        // NCP z-precision on the active block:
        //   P_z = I + diag(sig) G_m diag(sig)  =  I + (sig sig') .* G_m
        // Reuses G_m instead of a second weighted crossprod.
        sig_red = sig_m(active_re);
        P_z_red = arma::eye(active_idx.size(), active_idx.size()) +
                  (sig_red * sig_red.t()) % G_m.submat(active_re, active_re);
        // Pb_z = diag(sig) X_active' r_z = sig .* (X_active' r_z)
        Pb_z_red = sig_red % (Xm.cols(active_re).t() * r_z);

        z_red = chol_sample_precision_cpp(P_z_red, Pb_z_red);
        z_m(active_re) = z_red;
      }

      z_c_out.slice(m).col(ip) = z_m;

      // beta deviation in original coords: bdev = sig .* z  (zero off active),
      // so Xm_sc * z == Xm * bdev — no scaled design copy needed.
      bdev = sig_m % z_m;
      r_mu = kappa_w(idx, arma::uvec({(unsigned int)j})) +
             om_m % c_j(idx) - om_m % (Xm * bdev);
      if (use_bart)
        r_mu -= om_m % f_bart(idx, arma::uvec({(unsigned int)ip}));

      P_mu_acc += G_m;
      Pb_mu_acc += Xmt * r_mu;
    }

    arma::mat P_mu = prior_P + P_mu_acc;
    P_mu.diag() += hs_prec_mu.col(ip);
    arma::vec Pb_mu = prior_Pb.col(ip) + Pb_mu_acc;

    mu_new.col(ip) = chol_sample_precision_cpp(P_mu, Pb_mu);

    // Recover beta_c from centered z and updated mu
    for (int m = 0; m < n_groups; m++) {
      arma::vec z_m = z_c_out.slice(m).col(ip);
      for (int v = 0; v < k; v++) {
        bool active_v = is_re(v) && (re_mask(v, m) > 0.0) && (is_intercept(v) == 1 || y_mask(ip, m) > 0.0);
        beta_c(v, ip, m) = active_v
                               ? mu_new(v, ip) + sigma_mat(v, ip) * re_support(v, m) * z_m(v)
                               : mu_new(v, ip);
      }
    }

    // ---- RE-ASIS interweave: redraw each per-(covariate,category-ip) RE scale in the
    // non-centered parameterization (dev = sigma * re_support * z_c, z_c fixed) via an MH step.
    // Conditional on the half-Cauchy auxiliary a_re, the scale prior is IG(1/2, 1/a_re), so the
    // PG-weighted likelihood proposal cancels -> IG-prior ratio + sigma^2->sigma Jacobian. The
    // rescaled beta_c flow into the subsequent CP precision draw (in R) -> the interweave.
    // Sequential per-covariate (RFip updated after each accept), mirroring the CLR/LNM step.
    // Guarded against extreme rescales (PG-overflow). Default off; on==off in mean at re_idx=1.
    if (re_asis) {
      arma::vec RFip(n, arma::fill::zeros), resid_mu(n, arma::fill::zeros);
      for (int m = 0; m < n_groups; m++) {
        IntegerVector iv = idx_vec[m];
        for (int t = 0; t < iv.size(); t++) {
          int i = iv[t] - 1; double rf = 0.0;
          for (int vv = 0; vv < k; vv++) if (is_re(vv)) rf += X(i, vv) * (beta_c(vv, ip, m) - mu_new(vv, ip));
          RFip[i] = rf;
          resid_mu[i] = kappa_w(i, j) / om_p[i] + c_j[i] - arma::dot(X.row(i), mu_new.col(ip));
        }
      }
      for (int v = 0; v < k; v++) {
        if (!is_re(v)) continue;
        double s_old = sigma_mat(v, ip); if (!std::isfinite(s_old) || s_old < 1e-8) continue;
        double A = 0.0, C = 0.0;
        for (int m = 0; m < n_groups; m++) {
          double sc = re_support(v, m) * z_c_out(v, ip, m); if (sc == 0.0) continue;
          IntegerVector iv = idx_vec[m];
          for (int t = 0; t < iv.size(); t++) { int i = iv[t] - 1; double dv = X(i, v) * sc;
            A += om_p[i] * dv * dv; C += om_p[i] * dv * (resid_mu[i] - (RFip[i] - s_old * dv)); }
        }
        if (!std::isfinite(A) || A <= 0) continue;
        double s_prop = R::rnorm(C / A, 1.0 / std::sqrt(A));
        if (!std::isfinite(s_prop) || s_prop <= 0.1 * s_old || s_prop >= 10.0 * s_old) continue;
        double b_ig = 1.0 / std::max(a_re(v, ip), 1e-12);
        double la = (-1.5 * std::log(s_prop * s_prop) - b_ig / (s_prop * s_prop) + std::log(s_prop))
                  - (-1.5 * std::log(s_old * s_old)   - b_ig / (s_old * s_old)   + std::log(s_old));
        if (std::isfinite(la) && std::log(R::runif(0, 1)) < la) {
          double f = s_prop / s_old;
          for (int m = 0; m < n_groups; m++) {
            double sc = re_support(v, m) * z_c_out(v, ip, m);
            beta_c(v, ip, m) = mu_new(v, ip) + s_prop * sc;
            IntegerVector iv = idx_vec[m];
            for (int t = 0; t < iv.size(); t++) { int i = iv[t] - 1; RFip[i] += (f - 1.0) * s_old * X(i, v) * sc; }
          }
        }
      }
    }
  }

  // Finalize return
  NumericVector beta_c_v(k * p * n_groups);
  std::copy(beta_c.begin(), beta_c.end(), beta_c_v.begin());
  beta_c_v.attr("dim") = Dimension(k, p, n_groups);

  NumericMatrix mu_out(k, p);
  std::copy(mu_new.begin(), mu_new.end(), mu_out.begin());

  NumericVector z_c_v(k * p * n_groups);
  std::copy(z_c_out.begin(), z_c_out.end(), z_c_v.begin());
  z_c_v.attr("dim") = Dimension(k, p, n_groups);

  return List::create(Named("beta_c") = beta_c_v, Named("mu") = mu_out,
                      Named("z_c") = z_c_v);
}



// ==========================================================================
// ADD THIS FUNCTION to mnlogit_gibbs_core.cpp
// (replaces / supplements update_re_precision for the half-Cauchy path)
//
// Half-Cauchy prior on σ = 1/√τ via the Wand et al. (2011) auxiliary
// variable representation:
//
//   τ | a  ~  Gamma(1/2,  a)          [shape, rate]
//   a       ~  Gamma(1/2,  1/A²)
//
// which implies marginally  σ ~ half-Cauchy(A).
//
// Gibbs updates (both conjugate, no MH needed):
//   a | τ^{t-1}  ~  Exponential( τ^{t-1} + 1/A² )
//                 ≡  Gamma( 1, τ^{t-1} + 1/A² )
//   τ | a^t, data ~  Gamma( 1/2 + N_active/2,  a^t + SS/2 )
//
// Arguments:
//   beta_c        k × p × n_groups cube  (group-level effects)
//   mu_pooled     k × p matrix           (pooled means)
//   re_idx        0-based index vector of random-effect rows
//   n_groups      number of groups
//   re_mask       k × n_groups  (structural activity mask)
//   y_mask        p × n_groups  (category-group activity mask)
//   is_intercept  k-vector of {0,1}
//   prec_prev     k × p  — precision from the PREVIOUS iteration  (τ^{t-1})
//   a_aux_prev    k × p  — auxiliary variable from previous iter   (a^{t-1})
//   re_scale_A    half-Cauchy scale A  (default 1.0; try 2.5 for diffuse)
// ==========================================================================
// [[Rcpp::export]]
List update_re_precision_hc(const arma::cube   &beta_c,
                            const arma::mat    &mu_pooled,
                            const arma::uvec   &re_idx,
                            int                 n_groups,
                            const arma::cube   &re_mask,
                            const arma::mat    &y_mask,
                            const arma::uvec   &is_intercept,
                            const arma::mat    &prec_prev,
                            const arma::mat    &a_aux_prev,
                            double              re_scale_A = 1.0) {

  int k = beta_c.n_rows;
  int p = beta_c.n_cols;

  arma::mat prec_out  = prec_prev;    // carry forward for inactive cells
  arma::mat sigma_out = 1.0 / arma::sqrt(prec_out + 1e-12);  // initialize consistently from prec_prev
  arma::mat a_aux_out = a_aux_prev;   // carry forward for inactive cells

  const double inv_A2 = 1.0 / (re_scale_A * re_scale_A);

  for (int ip = 0; ip < p; ip++) {
    for (arma::uword ri = 0; ri < re_idx.n_elem; ri++) {
      int v = static_cast<int>(re_idx(ri));

      // ---- Accumulate sufficient statistics --------------------------------
      double ss       = 0.0;
      double n_active = 0.0;

      for (int m = 0; m < n_groups; m++) {
        if (re_mask(v, ip, m) > 0.5 &&
            (is_intercept(v) == 1 || y_mask(ip, m) > 0.5)) {
          double diff = beta_c(v, ip, m) - mu_pooled(v, ip);
          ss       += diff * diff;
          n_active += 1.0;
        }
      }

      // ---- Inactive cell: keep previous values, do not update -------------
      if (n_active < 1.0) {
        sigma_out(v, ip) = 1.0 / std::sqrt(prec_prev(v, ip) + 1e-12);
        continue;
      }

      // ---- Step 1: update auxiliary  a | τ^{t-1} --------------------------
      //   a ~ Gamma(1,  τ^{t-1} + 1/A²)  ≡  Exponential(τ^{t-1} + 1/A²)
      double tau_prev = std::max(prec_prev(v, ip), 1e-8);
      double rate_a   = tau_prev + inv_A2;
      double a_new    = R::rgamma(1.0, 1.0 / rate_a);
      a_new           = std::max(a_new, 1e-12);
      a_aux_out(v, ip) = a_new;

      // ---- Step 2: update precision  τ | a^t, data  -----------------------
      //   τ ~ Gamma( 0.5 + N_active/2,  a^t + SS/2 )
      double shape_tau = 0.5 + n_active / 2.0;
      double rate_tau  = a_new + ss / 2.0;
      if (!std::isfinite(rate_tau) || rate_tau <= 0.0) rate_tau = inv_A2;

      double prec = R::rgamma(shape_tau, 1.0 / rate_tau);
      prec = std::min(std::max(prec, 1e-4), 1e6);

      prec_out(v, ip)  = prec;
      sigma_out(v, ip) = 1.0 / std::sqrt(prec + 1e-12);
    }
  }

  return List::create(Named("prec")  = prec_out,
                      Named("sigma") = sigma_out,
                      Named("a_aux") = a_aux_out);
}


// ==========================================================================
// NEW FUNCTION: Symmetric pooled RE precision update
// ==========================================================================
// [[Rcpp::export]]
List update_re_precision_hc_sym(const arma::cube   &beta_c,
                                const arma::mat    &mu_pooled,
                                const arma::uvec   &re_idx,
                                int                 n_groups,
                                const arma::cube   &re_mask,
                                const arma::mat    &y_mask,
                                const arma::uvec   &is_intercept,
                                const arma::mat    &prec_prev,
                                const arma::mat    &a_aux_prev,
                                int                 p_all,
                                double              re_scale_A = 1.0,
                                Rcpp::Nullable<Rcpp::IntegerVector> block_id_opt = R_NilValue,
                                Rcpp::Nullable<Rcpp::IntegerVector> block_size_opt = R_NilValue,
                                bool                re_regularize = false,
                                double              slab_c2 = 100.0,
                                Rcpp::Nullable<Rcpp::NumericMatrix> re_support_opt = R_NilValue) {

  int k = beta_c.n_rows;
  int p = beta_c.n_cols; // p = p_all - 1

  // Support prior consistency: the RE DRAW (gibbs_step_re_ncp) scales the prior SD by
  // re_support_mat (per-group participation factor rs<=1), so d_{.,m} ~ N(0,(sigma*rs_m)^2).
  // To infer the BASE sigma consistently, the variance ss must be whitened: weight each group's
  // ss by 1/rs_m^2 (since (d/rs) ~ N(0,sigma^2)). Without this the variance step under-accounts
  // for the support shrinkage (the pre-existing inconsistency). rs default 1 -> no change.
  bool have_support = re_support_opt.isNotNull();
  arma::mat re_support;
  if (have_support) re_support = Rcpp::as<arma::mat>(re_support_opt);  // k x n_groups

  // Regularised horseshoe via a CONSISTENT slice (Finnish horseshoe done properly).
  // Effective variance v~ = c2*s2/(c2+s2) (s2=1/tau_raw) <=> effective precision
  // tau_eff = tau_raw + 1/c2. The REs use tau_eff EVERYWHERE; we infer tau_raw (the
  // half-Cauchy variable), whose conditional is NOT Gamma -> 1-D slice on log tau_raw.
  // We STORE tau_eff in prec_out (so the RE draw / sigma_mat use the regularised precision)
  // and recover tau_raw = tau_eff - 1/c2 for the aux. The collapsed path uses the SAME
  // tau_eff in its marginal+refill, so both target one well-behaved posterior (exp9). When
  // off (inv_c2=0) this is the exact unregularised Gamma -> identical to before.
  const double inv_A2 = 1.0 / (re_scale_A * re_scale_A);
  const double inv_c2 = (re_regularize && slab_c2 > 0.0 && std::isfinite(slab_c2)) ? 1.0 / slab_c2 : 0.0;
  // draw the (effective) RE precision given ss, n_active (na), aux a_new, and previous tau_eff
  auto draw_tau_eff = [&](double ss, double na, double a_new, double tau_prev_eff) -> double {
    double pr;
    if (inv_c2 <= 0.0) {                                   // exact unregularised Gamma conditional
      double rate = a_new + 0.5 * ss;
      if (!std::isfinite(rate) || rate <= 0.0) rate = inv_A2;
      pr = R::rgamma(0.5 + 0.5 * na, 1.0 / rate);
    } else {                                               // regularised: slice log(tau_raw), store tau_eff
      auto lp = [&](double lt) { double tr = std::exp(lt), te = tr + inv_c2;
        return 0.5 * na * std::log(te) - 0.5 * te * ss - 0.5 * lt - a_new * tr + lt; };
      double lt = std::log(std::max(tau_prev_eff - inv_c2, 1e-8));
      double y0 = lp(lt) - R::exp_rand();
      double Lb = lt - R::unif_rand(), Rb = Lb + 1.0; int gi = 0;
      while (lp(Lb) > y0 && gi < 80) { Lb -= 1.0; ++gi; } gi = 0;
      while (lp(Rb) > y0 && gi < 80) { Rb += 1.0; ++gi; }
      double ln = lt;
      for (int s2i = 0; s2i < 100; ++s2i) { double q = Lb + R::unif_rand() * (Rb - Lb);
        if (lp(q) > y0) { ln = q; break; }
        if (q < lt) Lb = q; else Rb = q; }
      pr = std::exp(ln) + inv_c2;
    }
    return std::min(std::max(pr, 1e-4), 1e6);
  };

  arma::mat prec_out  = prec_prev;
  arma::mat sigma_out = 1.0 / arma::sqrt(prec_out + 1e-12);
  arma::mat a_aux_out = a_aux_prev;

  std::vector<int> blk_id(k, -1);
  std::vector<int> blk_size;
  if (block_id_opt.isNotNull() && block_size_opt.isNotNull()) {
    Rcpp::IntegerVector b_id(block_id_opt);
    Rcpp::IntegerVector b_sz(block_size_opt);
    for (int i = 0; i < k; i++) blk_id[i] = b_id[i];
    for (int i = 0; i < b_sz.size(); i++) blk_size.push_back(b_sz[i]);
  }
  std::vector<bool> block_processed(blk_size.size(), false);

  for (arma::uword ri = 0; ri < re_idx.n_elem; ri++) {
    int v = static_cast<int>(re_idx(ri));

    int b_id = blk_id[v];
    if (b_id >= 0) {
      if (block_processed[b_id]) continue;
      block_processed[b_id] = true;
      
      double ss       = 0.0;
      double n_active = 0.0;
      int Kb = blk_size[b_id];
      
      std::vector<int> block_vars;
      for (arma::uword rj = 0; rj < re_idx.n_elem; rj++) {
        if (blk_id[re_idx(rj)] == b_id) {
          block_vars.push_back(re_idx(rj));
        }
      }
      
      for (int m = 0; m < n_groups; m++) {
        int num_active_categories = 0;
        std::vector<bool> ip_active(p, false);
        for (int ip = 0; ip < p; ip++) {
          bool active = false;
          for (int bv : block_vars) {
            if (re_mask(bv, ip, m) > 0.5 && (is_intercept(bv) == 1 || y_mask(ip, m) > 0.5)) {
              active = true;
              break;
            }
          }
          if (active) {
            ip_active[ip] = true;
            num_active_categories++;
          }
        }
        
        if (num_active_categories > 0) {
          std::vector<double> row_sums(block_vars.size(), 0.0);
          std::vector<double> col_sums(p_all, 0.0);
          double total_sum = 0.0;
          
          for (size_t iv = 0; iv < block_vars.size(); iv++) {
            int bv = block_vars[iv];
            for (int ip = 0; ip < p; ip++) {
              if (ip_active[ip]) {
                double diff = beta_c(bv, ip, m) - mu_pooled(bv, ip);
                row_sums[iv] += diff;
                col_sums[ip] += diff;
                total_sum += diff;
              }
            }
          }
          
          double mean_total = total_sum / (Kb * p_all);
          std::vector<double> mean_row(block_vars.size(), 0.0);
          for (size_t iv = 0; iv < block_vars.size(); iv++) mean_row[iv] = row_sums[iv] / p_all;
          double mean_row_dropped = 0.0 / p_all; 
          
          std::vector<double> mean_col(p_all, 0.0);
          for (int ip = 0; ip < p_all; ip++) mean_col[ip] = col_sums[ip] / Kb;
          
          for (size_t iv = 0; iv < block_vars.size(); iv++) {
            int bv = block_vars[iv];
            for (int ip = 0; ip < p; ip++) {
              if (ip_active[ip]) {
                double diff = beta_c(bv, ip, m) - mu_pooled(bv, ip);
                double d_centered = diff - mean_row[iv] - mean_col[ip] + mean_total;
                ss += d_centered * d_centered;
              }
            }
            double diff = 0.0;
            double d_centered = diff - mean_row[iv] - mean_col[p_all - 1] + mean_total;
            ss += d_centered * d_centered;
          }
          
          for (int ip = 0; ip < p; ip++) {
            if (ip_active[ip]) {
              double diff = 0.0;
              double d_centered = diff - mean_row_dropped - mean_col[ip] + mean_total;
              ss += d_centered * d_centered;
            }
          }
          double diff_d = 0.0;
          double d_centered_d = diff_d - mean_row_dropped - mean_col[p_all - 1] + mean_total;
          ss += d_centered_d * d_centered_d;
          
          n_active += num_active_categories * (Kb - 1);
        }
      }
      
      if (n_active < 1.0) {
        for (int bv : block_vars) {
          for (int ip = 0; ip < p; ip++) {
            sigma_out(bv, ip) = 1.0 / std::sqrt(prec_prev(bv, ip) + 1e-12);
          }
        }
        continue;
      }
      
      double tau_prev_sum = 0.0;
      for (int bv : block_vars) {
        for (int ip = 0; ip < p; ip++) {
          tau_prev_sum += prec_prev(bv, ip);
        }
      }
      double tau_prev = std::max(tau_prev_sum / (block_vars.size() * p), 1e-8);
      double rate_a   = std::max(tau_prev - inv_c2, 1e-8) + inv_A2;   // aux is on tau_raw = tau_eff - 1/c2
      double a_new    = std::max(R::rgamma(1.0, 1.0 / rate_a), 1e-12);
      double prec     = draw_tau_eff(ss, n_active, a_new, tau_prev);  // Gamma if off, slice if regularised

      for (int bv : block_vars) {
        for (int ip = 0; ip < p; ip++) {
          prec_out(bv, ip)  = prec;
          sigma_out(bv, ip) = 1.0 / std::sqrt(prec + 1e-12);
          a_aux_out(bv, ip) = a_new;
        }
      }
    } else {
      double ss       = 0.0;
      double n_active = 0.0;

      for (int m = 0; m < n_groups; m++) {
        double sum_diff = 0.0;
        int num_active_in_group = 0;
        
        for (int ip = 0; ip < p; ip++) {
          if (re_mask(v, ip, m) > 0.5 && (is_intercept(v) == 1 || y_mask(ip, m) > 0.5)) {
            double diff = beta_c(v, ip, m) - mu_pooled(v, ip);
            sum_diff += diff;
            num_active_in_group++;
          }
        }

        if (num_active_in_group > 0) {
          double mean_diff = sum_diff / p_all;
          double ss_m = 0.0;
          for (int ip = 0; ip < p; ip++) {
            if (re_mask(v, ip, m) > 0.5 && (is_intercept(v) == 1 || y_mask(ip, m) > 0.5)) {
              double diff = beta_c(v, ip, m) - mu_pooled(v, ip);
              double d_centered = diff - mean_diff;
              ss_m += d_centered * d_centered;
            }
          }
          double d_baseline = 0.0 - mean_diff;
          ss_m += d_baseline * d_baseline;

          double inv_rs2 = 1.0;                                   // support-consistent whitening
          if (have_support) { double rsv = re_support(v, m); if (rsv > 1e-12) inv_rs2 = 1.0 / (rsv * rsv); }
          ss += ss_m * inv_rs2;

          n_active += num_active_in_group;
        }
      }

      if (n_active < 1.0) {
        for (int ip = 0; ip < p; ip++) {
          sigma_out(v, ip) = 1.0 / std::sqrt(prec_prev(v, ip) + 1e-12);
        }
        continue;
      }

      double tau_prev_sum = 0.0;
      for (int ip = 0; ip < p; ip++) {
        tau_prev_sum += prec_prev(v, ip);
      }
      double tau_prev = std::max(tau_prev_sum / p, 1e-8);
      double rate_a   = std::max(tau_prev - inv_c2, 1e-8) + inv_A2;   // aux is on tau_raw = tau_eff - 1/c2
      double a_new    = std::max(R::rgamma(1.0, 1.0 / rate_a), 1e-12);
      double prec     = draw_tau_eff(ss, n_active, a_new, tau_prev);  // Gamma if off, slice if regularised

      for (int ip = 0; ip < p; ip++) {
        prec_out(v, ip)  = prec;
        sigma_out(v, ip) = 1.0 / std::sqrt(prec + 1e-12);
        a_aux_out(v, ip) = a_new;
      }
    }
  }

  return List::create(Named("prec")  = prec_out,
                      Named("sigma") = sigma_out,
                      Named("a_aux") = a_aux_out);
}


