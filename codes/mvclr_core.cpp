// mvclr_core.cpp — RcppArmadillo accelerators for the multivariate CLR Gibbs.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Matrix-normal draw:  B ~ MN(M, U, Sigma)  with  P = U^{-1}  (P = prior_prec + data),
//   M = P^{-1} RHS,  draw = M + L^{-T} Z Schol,   P = L L',  Schol' Schol = Sigma.
// Returns a (d x m) draw.  RHS is (d x m), Schol is the upper-chol of Sigma (m x m).
// [[Rcpp::export]]
arma::mat mn_draw_cpp(const arma::mat& P, const arma::mat& RHS, const arma::mat& Schol) {
  arma::mat L = arma::chol(P, "lower");                       // P = L L'
  arma::mat M = arma::solve(arma::trimatu(L.t()),
                            arma::solve(arma::trimatl(L), RHS)); // P^{-1} RHS
  arma::mat Z = arma::randn<arma::mat>(P.n_rows, RHS.n_cols);
  arma::mat LinvtZ = arma::solve(arma::trimatu(L.t()), Z);    // L^{-T} Z  (left factor of U)
  return M + LinvtZ * Schol;
}

// Per-group random-effect block draw (random intercept OR slopes).
//   For each group g: b_g ~ MN( (Xg'Xg + Prc)^{-1} Xg'Rg , (Xg'Xg+Prc)^{-1}, Sigma ).
//   Xre  : (n x q) random-effect design (subset of columns), Resid : (n x m) = W - X*B_fixed,
//   gstart/gsize : 0-based row ranges per group (rows assumed grouped/contiguous),
//   Prc  : (q x q) RE prior precision (diag),  Schol : upper-chol of Sigma.
// Returns (G x q x m) flattened as a (G*q) x m matrix (row block g = rows g*q..g*q+q-1).
// [[Rcpp::export]]
arma::mat re_block_draw_cpp(const arma::mat& Xre, const arma::mat& Resid,
                            const arma::uvec& gstart, const arma::uvec& gsize,
                            const arma::mat& Prc, const arma::mat& Schol) {
  unsigned int G = gstart.n_elem, q = Xre.n_cols, m = Resid.n_cols;
  arma::mat out(G * q, m, arma::fill::zeros);
  for (unsigned int g = 0; g < G; ++g) {
    arma::uvec rows = arma::regspace<arma::uvec>(gstart[g], gstart[g] + gsize[g] - 1);
    arma::mat Xg = Xre.rows(rows);
    arma::mat P  = Xg.t() * Xg + Prc;
    arma::mat RHS = Xg.t() * Resid.rows(rows);
    out.rows(g * q, g * q + q - 1) = mn_draw_cpp(P, RHS, Schol);
  }
  return out;
}
