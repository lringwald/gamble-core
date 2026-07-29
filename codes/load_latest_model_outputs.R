# =========================================================================
# Helper Script: Dynamically Load Latest MNL Model Outputs
# =========================================================================

#' Find the most recently modified file matching a specific pattern
#' @param dir_path Directory to search in
#' @param pattern Regex pattern to match files
#' @return Absolute path of the newest file, or NULL if none found
get_latest_file <- function(dir_path, pattern) {
  files <- list.files(path = dir_path, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) return(NULL)
  
  # Sort files by modification time (newest first)
  info <- file.info(files)
  latest_file <- rownames(info[order(info$mtime, decreasing = TRUE), ])[1]
  
  return(latest_file)
}

#' Load the latest covariate matrix and posterior beta means for a resolution
#' @param downscale_reso Resolution string (e.g., "5km", "20km")
#' @param base_path Optional base path to the high_res_prior_models directory
#' @return A list containing the loaded dat_pixel dataframe and beta_means dataframe
load_latest_model_outputs <- function(downscale_reso = "20km", base_path = NULL) {
  
  if (is.null(base_path)) {
    # Default path setup
    user_home <- Sys.getenv("HOME")
    base_path <- file.path(user_home, "Library/CloudStorage/OneDrive-IIASA/R/LAMASUS_high_res_prior_models")
  }
  
  output_dir <- file.path(base_path, "output")
  
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist: ", output_dir)
  }
  
  # Define regex patterns for the files
  # Matches e.g., dat_pixel_5km_2026-05-08.rds (removes ^ to allow matching in subfolders)
  xmat_pattern <- paste0("dat_pixel_", downscale_reso, "_[0-9]{4}-[0-9]{2}-[0-9]{2}\\.rds$")
  beta_pattern <- paste0("mnl_pixel_beta_means_", downscale_reso, "_[0-9]{4}-[0-9]{2}-[0-9]{2}\\.csv$")
  
  # Grab the most recent files
  xmat_path <- get_latest_file(output_dir, xmat_pattern)
  beta_path <- get_latest_file(output_dir, beta_pattern)
  
  # Validation
  if (is.null(xmat_path)) stop(sprintf("Could not find dat_pixel file for %s in %s", downscale_reso, output_dir))
  if (is.null(beta_path)) stop(sprintf("Could not find beta means file for %s in %s", downscale_reso, output_dir))
  
  cat(sprintf("Loading latest covariates from: %s\n", basename(xmat_path)))
  dat_pixel <- readRDS(xmat_path)
  
  cat(sprintf("Loading latest beta means from: %s\n", basename(beta_path)))
  beta_means <- read.csv(beta_path, stringsAsFactors = FALSE)
  
  return(list(
    dat_pixel = dat_pixel,
    beta_means = beta_means,
    xmat_path = xmat_path,
    beta_path = beta_path
  ))
}

# Example Usage:
# source("codes/load_latest_model_outputs.R")
# latest_data <- load_latest_model_outputs(downscale_reso = "5km")
# dat_pixel <- latest_data$dat_pixel
# beta_means <- latest_data$beta_means
