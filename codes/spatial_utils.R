#' Efficient Focal Statistics on Data Frames with Distance Decay
#'
#' @description Calculates neighborhood statistics for spatial grids stored as data frames.
#' Supports distance-based weighting (decay).
#'
#' @param df A data frame containing spatial grid data.
#' @param x_id_col Name of the column representing the x-coordinate or x-index.
#' @param y_id_col Name of the column representing the y-coordinate or y-index.
#' @param value_cols Vector of names of the columns to calculate statistics for.
#' @param w Neighborhood window size (must be odd, e.g., 3, 5) or a matrix of weights.
#' @param fun The function to apply (e.g., "mean", "sum", "var"). 
#' @param na.rm Logical. Should missing values be removed?
#' @param fill Number. Value to fill for cells outside the grid (default is NA).
#' @param include_center Logical. Should the center cell be included? Default is FALSE.
#' @param decay_type String. Type of distance decay: "none", "inverse", or "exponential".
#' @param decay_param Numeric. Parameter for decay (power for "inverse", rate for "exponential").
#'
#' @return A data frame with the original coordinates and the calculated statistics.
#' @export
focal_df <- function(df, x_id_col, y_id_col, value_cols, w = 3, fun = "mean", 
                     na.rm = TRUE, fill = NA, include_center = FALSE,
                     decay_type = "none", decay_param = 1) {
  require(terra)
  require(dplyr)
  original_fun <- fun

  # 1. Prepare the weight matrix
  if (is.numeric(w) && length(w) == 1) {
    if (w %% 2 == 0) stop("Window size 'w' must be an odd integer.")
    
    # Coordinates of cells relative to center
    c_idx <- (w + 1) / 2
    row_idx <- rep(1:w, each = w)
    col_idx <- rep(1:w, times = w)
    
    # Euclidean distance from center
    dist_mat <- matrix(sqrt((row_idx - c_idx)^2 + (col_idx - c_idx)^2), nrow = w, ncol = w)
    
    if (decay_type == "none") {
      w_mat <- matrix(1, w, w)
    } else if (decay_type == "inverse") {
      # Use a small epsilon to avoid division by zero if center is included
      eps <- 1e-6
      w_mat <- 1 / (dist_mat + eps)^decay_param
    } else if (decay_type == "exponential") {
      w_mat <- exp(-decay_param * dist_mat)
    } else {
      stop("Invalid decay_type. Choose 'none', 'inverse', or 'exponential'.")
    }

    if (!include_center) {
      w_mat[c_idx, c_idx] <- 0
    }
    
    # If fun is "mean", normalize w_mat so it sums to 1 and use "sum" focal operation
    # This is more robust than relying on terra's internal "mean" with a weight matrix
    if (fun == "mean") {
      w_sum <- sum(w_mat, na.rm = TRUE)
      if (w_sum > 0) w_mat <- w_mat / w_sum
      fun <- "sum"
    }
  } else if (is.matrix(w)) {
    w_mat <- w
  } else {
    stop("w must be a single numeric value (window size) or a matrix.")
  }

  # 2. Convert data frame to terra::SpatRaster
  r_cols <- c(x_id_col, y_id_col, value_cols)
  r <- terra::rast(df[, r_cols], type = "xyz")

  # 3. Apply focal operation
  if (is.character(fun) && !fun %in% c("sum", "mean", "min", "max", "modal", "median")) {
    if (exists(fun, mode = "function")) {
      fun_obj <- get(fun, mode = "function")
    } else {
      stop(paste("Function", fun, "not found."))
    }
  } else {
    fun_obj <- fun
  }

  r_focal <- terra::focal(
    r,
    w = w_mat,
    fun = fun_obj,
    na.rm = na.rm,
    fillvalue = fill
  )

  # 4. Convert back to data frame
  df_out <- terra::as.data.frame(r_focal, xy = TRUE)
  
  # 5. Rename columns
  # Use original_fun for suffix even if we used sum internally for mean calculation
  suffix <- if (decay_type == "none") paste0("_focal_", original_fun) else paste0("_focal_", original_fun, "_", decay_type)
  colnames(df_out) <- c(x_id_col, y_id_col, paste0(value_cols, suffix))

  return(df_out)
}
