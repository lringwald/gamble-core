#!/usr/bin/env Rscript
# =============================================================================
# diagnose_posterior.R — per-block convergence for a streamed posterior directory
# =============================================================================
#   Rscript diagnose_posterior.R output/bart_gate/posterior_bart
#
# Reads every chain written by save_posterior_to_disk and reports split-Rhat and bulk ESS for
# each STRUCTURAL ELEMENT separately, because they fail differently: the location parameters
# (mu, b_g) routinely mix two orders of magnitude better than the variance/shrinkage block, and
# a single averaged number hides that.
#
# Needs >= 2 chains for Rhat. With one chain it reports ESS only and says so.
# =============================================================================
suppressMessages({library(Rcpp); library(RcppArmadillo); library(posterior)})
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript diagnose_posterior.R <posterior_dir>")
DIR <- args[1]
if (!dir.exists(DIR)) stop("not a directory: ", DIR)
source("codes/mnl_aux_func.R"); source("codes/mnlogit_rcpp_sym.R")

ids <- sort(unique(as.integer(sub(".*_chain_([0-9]+)\\.qs$", "\\1",
        list.files(DIR, pattern = "^posterior_batch_.*\\.qs$")))))
if (!length(ids)) stop("no posterior batches in ", DIR)
cat(sprintf("\n%s\n%s\n%s\n", strrep("=",78), sprintf("POSTERIOR DIAGNOSTICS  %s", DIR), strrep("=",78)))
cat(sprintf("chains found: %s\n", paste(ids, collapse=", ")))
ch <- lapply(ids, function(i) recover_mnlogit_posterior(DIR, chain_id = i))

cap <- function(x, n = 500L) { if (is.null(x)) return(NULL)
  if (ncol(x) <= n) return(x); set.seed(1); x[, sort(sample(ncol(x), n)), drop = FALSE] }
as_mat <- function(z) { if (is.null(z)) return(NULL)
  if (is.vector(z)) return(matrix(z, ncol = 1))
  d <- dim(z); if (length(d) == 2L) t(z) else matrix(aperm(z, c(length(d), seq_len(length(d)-1))), nrow = d[length(d)]) }

BLOCKS <- list(
  "log_lik"      = function(f) as_mat(f$post_log_lik),
  "mu (FE)"      = function(f) cap(as_mat(f$postb_pooled)),
  "b_g (RE)"     = function(f) cap(as_mat(f$postb_total)),
  "sigma_re"     = function(f) cap(as_mat(f$post_sigma_re)),
  "kappa (HS)"   = function(f) cap(as_mat(f$horseshoe$post_kappa_pooled)),
  "delta (alt)"  = function(f) as_mat(f$post_delta),
  "slab_c2 (RE)" = function(f) as_mat(f$post_slab_c2),
  "re_tau"       = function(f) as_mat(f$post_re_tau),
  "kappa_v"      = function(f) as_mat(f$post_kappa))

cat(sprintf("\n%-14s %6s | %-26s | %-18s\n", "block", "n", "Rhat  med / max / >1.05", "ESS_bulk med / min"))
cat(sprintf("%s\n", strrep("-", 78)))
res <- list()
for (nm in names(BLOCKS)) {
  L <- lapply(ch, BLOCKS[[nm]])
  if (any(vapply(L, is.null, TRUE))) { cat(sprintf("%-14s   -- not present in this fit\n", nm)); next }
  n <- min(vapply(L, nrow, 1L)); p <- ncol(L[[1]])
  if (n < 4L) { cat(sprintf("%-14s   -- only %d draws\n", nm, n)); next }
  A <- array(unlist(lapply(L, function(m) m[seq_len(n), , drop = FALSE])), c(n, p, length(L)))
  d <- posterior::as_draws_array(aperm(A, c(1, 3, 2)))
  s <- posterior::summarise_draws(d, "rhat", "ess_bulk")
  s <- s[is.finite(s$ess_bulk), ]
  if (!nrow(s)) { cat(sprintf("%-14s   -- degenerate (constant)\n", nm)); next }
  rh <- s$rhat[is.finite(s$rhat)]
  cat(sprintf("%-14s %6d | %5s %6s %7s | %8.0f %8.0f\n", nm, nrow(s),
      if (length(rh)) sprintf("%.3f", median(rh)) else "  -",
      if (length(rh)) sprintf("%.3f", max(rh)) else "  -",
      if (length(rh)) sprintf("%4.0f%%", 100*mean(rh > 1.05)) else "  -",
      median(s$ess_bulk), min(s$ess_bulk)))
  res[[nm]] <- s
}
if (length(ids) < 2L)
  cat("\nNOTE: one chain only -- Rhat is undefined and shown as '-'. Rerun with BG_CHAINS>=2 for it.\n")

# BART surface: not a sampled parameter, so summarise how much it MOVES per class instead
if (!is.null(ch[[1]]$tree_store)) {
  cat(sprintf("\n%s\nBART surface\n%s\n", strrep("-",78), strrep("-",78)))
  cat(sprintf("  tree draws per chain: %s\n", paste(vapply(ch, function(f) length(f$tree_store), 1L), collapse=", ")))
  cat(sprintf("  ensembles per draw:   %d %s\n", length(ch[[1]]$tree_store[[1]]),
      if (isTRUE(ch[[1]]$bart_symmetric)) "(symmetric: one per class)" else "(one per NON-baseline class)"))
  cat("  f is a deterministic function of the trees, not a sampled scalar -- for its per-class\n")
  cat("  spread see the sd(f) column in the gate output, or reconstruct_bart_f_mean() on new X.\n")
}
cat(sprintf("\nRule of thumb: Rhat > 1.05 or ESS < 100 on log_lik means the JOINT chain has not mixed,\n"))
cat("even when mu looks healthy -- per-parameter medians flatter a chain that has not converged.\n")
