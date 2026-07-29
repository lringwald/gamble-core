// lnm_core.cpp — RcppArmadillo accelerators for the logistic-normal multinomial.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Per-pixel overdispersion draw: u_i ~ N(P^{-1} b_i, P^{-1}),  P = diag(om_i) + SigInv.
//   om, rb : n x m   (rb = z - X B - RE residual);  SigInv : m x m.  Returns n x m.
// This is the 8.6 s/iter bottleneck in pure R (n tiny Cholesky decompositions).
// [[Rcpp::export]]
arma::mat lnm_u_draw_cpp(const arma::mat& om, const arma::mat& rb, const arma::mat& SigInv) {
  unsigned int n = om.n_rows, m = om.n_cols;
  arma::mat U(n, m);
  for (unsigned int i = 0; i < n; ++i) {
    arma::mat P = SigInv;
    P.diag() += om.row(i).t();
    arma::mat L = arma::chol(P, "lower");
    arma::vec bi = (om.row(i) % rb.row(i)).t();
    arma::vec mu = arma::solve(arma::trimatu(L.t()), arma::solve(arma::trimatl(L), bi));
    U.row(i) = (mu + arma::solve(arma::trimatu(L.t()), arma::randn<arma::vec>(m))).t();
  }
  return U;
}

// Country random effects (intercept OR slopes), per (group g, category j), with
// category-specific PG weights. For each g,j: b_gj ~ N over a q-dim weighted
// Gaussian regression of resid on Xre. Returns (G*q) x m (group g block = rows g*q..).
//   Xre   : n x q (group-sorted),  resid : n x m (= z - X B - u),  om : n x m,
//   gstart/gsize : 0-based group row ranges,  tau_re2 : RE prior variance (per category-col).
// [[Rcpp::export]]
arma::mat lnm_re_draw_cpp(const arma::mat& Xre, const arma::mat& resid, const arma::mat& om,
                          const arma::uvec& gstart, const arma::uvec& gsize, const arma::vec& tau_re2) {
  unsigned int G = gstart.n_elem, q = Xre.n_cols, m = resid.n_cols;
  arma::mat out(G * q, m, arma::fill::zeros);
  for (unsigned int g = 0; g < G; ++g) {
    unsigned int r0 = gstart[g], r1 = gstart[g] + gsize[g] - 1;   // contiguous (rows are group-sorted)
    arma::mat Xg = Xre.rows(r0, r1);
    for (unsigned int j = 0; j < m; ++j) {
      arma::vec wj = om.col(j).subvec(r0, r1);
      arma::mat P = Xg.t() * (Xg.each_col() % wj);
      P.diag() += 1.0 / tau_re2[j];
      arma::vec b = Xg.t() * (wj % resid.col(j).subvec(r0, r1));
      arma::mat L = arma::chol(P, "lower");
      arma::vec mu = arma::solve(arma::trimatu(L.t()), arma::solve(arma::trimatl(L), b));
      out.submat(g * q, j, g * q + q - 1, j) =
        mu + arma::solve(arma::trimatu(L.t()), arma::randn<arma::vec>(q));
    }
  }
  return out;
}

// Per-category fixed-effect draw: B_j ~ N(P^{-1} X' (om_j * r_j), P^{-1}),
//   P = X' diag(om_j) X + diag(prior_prec).  X : n x p,  r,om : n x m.  Returns p x m.
// [[Rcpp::export]]
arma::mat lnm_beta_draw_cpp(const arma::mat& X, const arma::mat& r, const arma::mat& om,
                            const arma::vec& prior_prec) {
  unsigned int p = X.n_cols, m = r.n_cols;
  arma::mat B(p, m);
  for (unsigned int j = 0; j < m; ++j) {
    arma::vec wj = om.col(j);
    arma::mat P = X.t() * (X.each_col() % wj);
    P.diag() += prior_prec;
    arma::vec b = X.t() * (wj % r.col(j));
    arma::mat L = arma::chol(P, "lower");
    arma::vec mu = arma::solve(arma::trimatu(L.t()), arma::solve(arma::trimatl(L), b));
    B.col(j) = mu + arma::solve(arma::trimatu(L.t()), arma::randn<arma::vec>(p));
  }
  return B;
}
