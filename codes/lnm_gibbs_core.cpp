// lnm_gibbs_core.cpp — FULL logistic-normal multinomial Gibbs loop in C++.
// PG via the large-h normal approximation (exact-accurate for levels N>>1; the
// LNM only runs on levels). Validates against the exact-PG R version.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

inline double pg_mean_c(double b, double c) {
  double ac = std::fabs(c); if (ac < 1e-6) return b * 0.25;
  return (b / (2.0 * c)) * std::tanh(c / 2.0);
}
inline double pg_var_c(double b, double c) {
  double ac = std::fabs(c); if (ac < 1e-6) return b / 24.0;
  double cc = std::min(ac, 30.0) * (c < 0 ? -1.0 : 1.0);          // guard overflow
  double ch = std::cosh(cc / 2.0);
  return std::max(b * (std::sinh(cc) - cc) / (4.0 * cc * cc * cc * ch * ch), 1e-12);
}

// [[Rcpp::export]]
List lnm_gibbs_core(const arma::mat& X, const arma::mat& Yb, const arma::vec& N,
                    const arma::mat& Xre, const arma::uvec& gstart, const arma::uvec& gsize,
                    const arma::uvec& hs_idx0, int niter, int nburn,
                    double A0, double slab_s2, double iw_df, const arma::mat& iw_scale,
                    double re_a, double re_b, bool use_re, bool use_hs,
                    int thin, unsigned int seed, bool use_asis, bool use_tempering = false, double tempering_T0 = 0.05) {
  arma::arma_rng::set_seed(seed);
  unsigned int n = X.n_rows, p = X.n_cols, m = Yb.n_cols, G = gstart.n_elem, q = Xre.n_cols;
  arma::mat kappa(n, m);
  for (unsigned int j = 0; j < m; ++j) kappa.col(j) = Yb.col(j) - N / 2.0;

  arma::mat B(p, m, arma::fill::zeros), u(n, m, arma::fill::zeros);
  arma::mat bre(std::max(G * q, 1u), m, arma::fill::zeros);
  arma::mat Sigma = arma::eye(m, m), SigInv = arma::eye(m, m);
  arma::vec tau_re2(m, arma::fill::ones), lambda2(p, arma::fill::ones), nuhs(p, arma::fill::ones);
  arma::vec prior_prec(p); prior_prec.fill(1.0 / (A0 * A0));
  double tau2 = 1, xi = 1, c2 = slab_s2;
  arma::mat postB_sum(p, m, arma::fill::zeros), postSig_sum(m, m, arma::fill::zeros);
  arma::mat bre_sum(std::max(G * q, 1u), m, arma::fill::zeros);
  arma::mat tau_re_sum(m, 1, arma::fill::zeros);
  arma::vec maxB(niter), od(niter);
  int n_keep = (niter - nburn) / thin; if (n_keep < 1) n_keep = 1;
  arma::mat Bdraws(p * m, n_keep, arma::fill::zeros);      // thinned, symmetric-able post draws
  arma::mat oddraws(m, n_keep, arma::fill::zeros);         // Sigma_od diagonal per kept draw
  int kept = 0;

  for (int it = 0; it < niter; ++it) {
    double temp_iter = 1.0;
    if (use_tempering && it <= nburn / 4 && nburn > 0) {
      temp_iter = tempering_T0 + (1.0 - tempering_T0) * ((double)it / (nburn / 4));
    }
    arma::vec N_iter = N * temp_iter;
    arma::mat kappa_iter = kappa * temp_iter;
    arma::mat RF(n, m, arma::fill::zeros);
    if (use_re) for (unsigned int g = 0; g < G; ++g) {
      unsigned int r0 = gstart[g], r1 = gstart[g] + gsize[g] - 1;
      RF.rows(r0, r1) = Xre.rows(r0, r1) * bre.rows(g * q, g * q + q - 1);
    }
    arma::mat eta = X * B + RF + u;
    arma::mat expe = arma::exp(arma::clamp(eta, -30.0, 30.0));
    arma::vec denom = 1.0 + arma::sum(expe, 1);

    arma::mat om(n, m), z(n, m);
    for (unsigned int j = 0; j < m; ++j) {
      arma::vec cj = arma::log(arma::clamp(denom - expe.col(j), 1e-12, arma::datum::inf));
      arma::vec psi = eta.col(j) - cj;
      for (unsigned int i = 0; i < n; ++i) {
        double o = R::rnorm(pg_mean_c(N_iter[i], psi[i]), std::sqrt(pg_var_c(N_iter[i], psi[i])));
        om(i, j) = std::max(o, 1e-9);
      }
      z.col(j) = kappa_iter.col(j) / om.col(j) + cj;
    }

    if (use_hs) for (arma::uword idx : hs_idx0) {
      double ve = c2 * tau2 * lambda2[idx] / (c2 + tau2 * lambda2[idx]);
      prior_prec[idx] = 1.0 / std::max(ve, 1e-12);
    }
    // B | .
    arma::mat rB = z - RF - u;
    for (unsigned int j = 0; j < m; ++j) {
      arma::vec wj = om.col(j); arma::mat P = X.t() * (X.each_col() % wj); P.diag() += prior_prec;
      arma::mat L = arma::chol(P, "lower"); arma::vec b = X.t() * (wj % rB.col(j));
      arma::vec mu = arma::solve(arma::trimatu(L.t()), arma::solve(arma::trimatl(L), b));
      B.col(j) = mu + arma::solve(arma::trimatu(L.t()), arma::randn<arma::vec>(p));
    }
    // RE | .
    if (use_re) {
      arma::mat rr = z - X * B - u;
      for (unsigned int g = 0; g < G; ++g) {
        unsigned int r0 = gstart[g], r1 = gstart[g] + gsize[g] - 1; arma::mat Xg = Xre.rows(r0, r1);
        for (unsigned int j = 0; j < m; ++j) {
          arma::vec wj = om.col(j).subvec(r0, r1); arma::mat P = Xg.t() * (Xg.each_col() % wj); P.diag() += 1.0 / tau_re2[j] + 1e-6;  // ridge: keep per-group RE system PD with many random slopes / small groups
          arma::vec b = Xg.t() * (wj % rr.col(j).subvec(r0, r1)); arma::mat L = arma::chol(P, "lower");
          arma::vec mu = arma::solve(arma::trimatu(L.t()), arma::solve(arma::trimatl(L), b));
          bre.submat(g * q, j, g * q + q - 1, j) = mu + arma::solve(arma::trimatu(L.t()), arma::randn<arma::vec>(q));
        }
      }
      for (unsigned int j = 0; j < m; ++j) { double ss = arma::dot(bre.col(j), bre.col(j));
        tau_re2[j] = 1.0 / R::rgamma(re_a + G * q / 2.0, 1.0 / (re_b + ss / 2.0)); }
      RF.zeros(); for (unsigned int g = 0; g < G; ++g) { unsigned int r0 = gstart[g], r1 = gstart[g] + gsize[g] - 1;
        RF.rows(r0, r1) = Xre.rows(r0, r1) * bre.rows(g * q, g * q + q - 1); }
    }
    // u | .
    arma::mat rb = z - X * B - RF;
    for (unsigned int i = 0; i < n; ++i) {
      arma::mat P = SigInv; P.diag() += om.row(i).t(); arma::mat L = arma::chol(P, "lower");
      arma::vec bi = (om.row(i) % rb.row(i)).t(); arma::vec mu = arma::solve(arma::trimatu(L.t()), arma::solve(arma::trimatl(L), bi));
      u.row(i) = (mu + arma::solve(arma::trimatu(L.t()), arma::randn<arma::vec>(m))).t();
    }
    // Sigma_od ~ IW
    Sigma = arma::iwishrnd(iw_scale + u.t() * u, iw_df + n); SigInv = arma::inv_sympd(Sigma);
    // FULL non-centered ASIS interweave: eps = L^{-1} u (ancillary), then redraw
    // Sigma via Cholesky-row regressions  r_.j = sum_{k<=j} L_jk eps_.k + noise(1/om_.j).
    // Gives Sigma m(m+1)/2 dof to move with eps fixed (scalar PX failed: 1 dof pinned by n).
    if (use_asis) {
      arma::mat Lc = arma::chol(Sigma, "lower");
      arma::mat eps = arma::solve(arma::trimatl(Lc), u.t()).t();          // n x m
      arma::mat rr = z - X * B - RF; arma::mat Ln(m, m, arma::fill::zeros);
      for (unsigned int j = 0; j < m; ++j) {
        arma::mat E = eps.cols(0, j); arma::vec wj = om.col(j);
        arma::mat P = E.t() * (E.each_col() % wj); P.diag() += 1e-6;       // weak prior (n large)
        arma::vec b = E.t() * (wj % rr.col(j)); arma::mat Lp = arma::chol(P, "lower");
        arma::vec mu = arma::solve(arma::trimatu(Lp.t()), arma::solve(arma::trimatl(Lp), b));
        arma::vec dr = mu + arma::solve(arma::trimatu(Lp.t()), arma::randn<arma::vec>(j + 1));
        if (dr[j] <= 1e-6) dr[j] = std::fabs(dr[j]) + 1e-6;                // positive diagonal -> PD
        Ln.submat(j, 0, j, j) = dr.t();
      }
      Sigma = Ln * Ln.t(); SigInv = arma::inv_sympd(Sigma);
      u = (Ln * eps.t()).t();                                             // consistent u for new Sigma
    }
    // horseshoe
    if (use_hs) {
      // shrinkage driven by the SYMMETRIC zero-sum norm (order-invariant; baseline-coded
      // norm would make shrinkage depend on which category is the reference)
      auto sym_ss = [&](arma::uword idx) { arma::rowvec br = B.row(idx);
        double mu_f = arma::accu(br) / (m + 1.0);                 // mean over K incl. baseline 0
        return arma::accu(arma::square(br - mu_f)) + mu_f * mu_f; };  // + (0 - mu_f)^2
      for (arma::uword idx : hs_idx0) { double ss = sym_ss(idx);
        nuhs[idx] = 1.0 / R::rgamma(1.0, 1.0 / (1.0 + 1.0 / lambda2[idx]));
        lambda2[idx] = std::min(std::max(1.0 / R::rgamma((m + 1) / 2.0, 1.0 / (1.0 / nuhs[idx] + ss / (2.0 * tau2))), 1e-10), 1e6); }
      double ssa = 0; for (arma::uword idx : hs_idx0) ssa += sym_ss(idx) / lambda2[idx];
      xi = 1.0 / R::rgamma(1.0, 1.0 / (1.0 + 1.0 / tau2));
      tau2 = std::min(std::max(1.0 / R::rgamma((hs_idx0.n_elem * m + 1) / 2.0, 1.0 / (1.0 / xi + ssa / 2.0)), 1e-10), 1e6);
    }
    maxB[it] = arma::abs(B).max(); od[it] = arma::mean(Sigma.diag());
    if (it >= nburn) {
      postB_sum += B; postSig_sum += Sigma;
      if (use_re) { bre_sum += bre; tau_re_sum += tau_re2; }
      if (((it - nburn) % thin == 0) && kept < n_keep) {
        Bdraws.col(kept) = arma::vectorise(B); oddraws.col(kept) = Sigma.diag(); ++kept;
      }
    }
  }
  return List::create(_["B_sum"] = postB_sum, _["Sig_sum"] = postSig_sum, _["nret"] = niter - nburn,
                      _["maxB"] = maxB, _["od"] = od, _["bre"] = bre,
                      _["bre_mean"] = bre_sum / (niter - nburn), _["tau_re2"] = tau_re_sum / (niter - nburn),
                      _["Bdraws"] = Bdraws.cols(0, kept - 1), _["oddraws"] = oddraws.cols(0, kept - 1));
}
