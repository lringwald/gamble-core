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

# ---- WHERE PACKAGES COME FROM ------------------------------------------------------------------
# RESPECT THE IMAGE'S OWN REPOSITORY. rocker/* point at Posit Package Manager and the r2u stack
# routes install.packages through apt via bspm; both install BINARIES in seconds. This file used to
# pass repos = "https://cloud.r-project.org/" explicitly, which overrides either of them and forces
# a SOURCE build of dbarts, RcppArmadillo, arrow and sf -- tens of minutes and a lot of memory, to
# arrive at the same packages.
#
# Only when nothing is configured do we choose: Posit Package Manager for the running distribution
# (binaries), and source CRAN as the last resort.
.repos <- getOption("repos")
.configured <- !is.null(.repos) && "CRAN" %in% names(.repos) &&
               nzchar(.repos[["CRAN"]]) && !identical(unname(.repos[["CRAN"]]), "@CRAN@")
if (!.configured) {
  .code <- tryCatch({
    os <- readLines("/etc/os-release", warn = FALSE)
    sub("^VERSION_CODENAME=", "", grep("^VERSION_CODENAME=", os, value = TRUE)[1])
  }, error = function(e) NA_character_)
  .code <- if (length(.code) && !is.na(.code)) gsub('"', "", .code) else NA_character_
  options(repos = if (!is.na(.code) && nzchar(.code) && .Platform$OS.type == "unix")
                    c(CRAN = sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest", .code))
                  else c(CRAN = "https://cloud.r-project.org/"))
}
# bspm asks for this in a container; without it the apt route falls back to source silently.
if (requireNamespace("bspm", quietly = TRUE)) options(bspm.sudo = TRUE)
message(">>> [init.R] repository: ", getOption("repos")[["CRAN"]],
        if (.configured) "  (from the image)" else "  (chosen here)")

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
  # No `repos` argument: use whatever was resolved above, so a binary repository stays in play.
  install.packages(missing, Ncpus = max(1L, parallel::detectCores()))
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
