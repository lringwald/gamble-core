# =============================================================================
# init.R — Predefined Stack Initialization for gamble-core
# =============================================================================
message(">>> [init.R] Checking and installing gamble-core dependencies...")

pkgs <- c(
  "Rcpp", "RcppArmadillo", "data.table", "qs2", "matrixStats",
  "posterior", "abind", "future", "future.apply", "progressr",
  "ggplot2", "scales", "base64enc", "arrow", "sf", "terra"
)

missing <- setdiff(pkgs, installed.packages()[, "Package"])
if (length(missing) > 0) {
  message(">>> Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org/", Ncpus = parallel::detectCores())
} else {
  message(">>> All required packages already installed.")
}
