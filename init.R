# =============================================================================
# init.R — the dependency contract for gamble-core
# =============================================================================
# ONE list, used by BOTH deployment paths, so they cannot drift apart:
#   * Accelerator / predefined R stack — runs this file directly
#   * Dockerfile                        — RUNs this file instead of its own list
#
# The set below is DERIVED FROM THE CODE, not curated by hand: every package a
# `library()` call loads at source time on the fit path. Those calls are hard, so
# a missing one is not a degraded run -- it is a crash at startup, before a single
# sweep. Re-derive with:
#   grep -rhoE "library\([A-Za-z][A-Za-z0-9.]*" codes drivers | sort -u
# =============================================================================
message(">>> [init.R] Checking and installing gamble-core dependencies...")

pkgs <- c(
  # --- sampler core: loaded by codes/mnlogit_rcpp_sym.R and codes/count_rcpp.R ---
  # These are the ones the previous list omitted. The sampler sources them at the
  # TOP of the file, so without them nothing runs at all.
  "Rcpp", "RcppArmadillo", "pg", "BayesLogit", "TruncatedNormal", "dbarts",
  "loo", "MASS", "Matrix", "matrixStats", "progressr", "qs2",
  # --- design assembly + drivers ---
  "data.table", "dplyr", "tidyr", "arrow", "sf", "terra",
  "future", "future.apply", "abind", "posterior",
  # --- reporting and diagnostics ---
  "ggplot2", "scales", "patchwork", "base64enc", "bayesplot", "spdep",
  # --- nested_cut store fingerprinting (guarded by requireNamespace, but the
  #     resumable store silently loses its config hash without it) ---
  "digest"
)

missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) {
  message(">>> Installing ", length(missing), " missing package(s): ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org/",
                   Ncpus = max(1L, parallel::detectCores()))
} else {
  message(">>> All required packages already installed.")
}

# VERIFY, don't assume. install.packages() signals a failed build as a WARNING and
# returns normally, so without this check a stack missing a system library (GDAL for
# sf, the Arrow C++ runtime for arrow) reports a clean init and then fails hours
# later inside a routine, with an error that points at the model instead of the image.
still <- setdiff(pkgs, rownames(installed.packages()))
if (length(still)) {
  stop("[init.R] these packages could not be installed: ", paste(still, collapse = ", "),
       "\n  Usually a missing SYSTEM library rather than a missing R package.",
       "\n  The Dockerfile apt-installs the ones this repo needs (arrow/parquet, BLAS/LAPACK);",
       "\n  a predefined stack has to supply them itself.", call. = FALSE)
}
message(">>> [init.R] ", length(pkgs), " package(s) present.")
