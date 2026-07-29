library(ggplot2)
library(dplyr)
library(tidyr)

#' MNL Parameter Heatplot
#' 
#' Creates a tiled heatplot where:
#' - Color represents the posterior median.
#' - Tile size represents the posterior credibility (distance from 50/50 sign flip).
#' 
#' @param draws A 3D array of posterior draws [k x p x draws] or a matrix/data frame
#' @param cov_names Character vector of covariate names
#' @param cat_names Character vector of category names
#' @param title Plot title
#' @param subtitle Plot subtitle
#' 
#' @return A ggplot object
MNL_parameter_heatplot <- function(draws, cov_names = NULL, cat_names = NULL, title = "Coefficient Heatplot", subtitle = "Tile size = Credibility, Color = median", exclude_intercept = FALSE, scaling_vec = NULL) {
  
  # Ensure draws is an array or handle summary data
  if (is.array(draws)) {
    # Extract names from dimnames if not provided
    if (is.null(cov_names)) cov_names <- dimnames(draws)[[1]]
    if (is.null(cat_names)) cat_names <- dimnames(draws)[[2]]
    
    # Fallback to indices if names are still missing
    if (is.null(cov_names)) cov_names <- paste0("V", 1:dim(draws)[1])
    if (is.null(cat_names)) cat_names <- paste0("C", 1:dim(draws)[2])

    # Compute summary statistics
    n_covs <- dim(draws)[1]
    n_cats <- dim(draws)[2]
    
    sum_vals <- apply(draws, 1:2, function(x) {
      p_val <- mean(x <= 0)
      c(median(x), abs(p_val - 0.5) * 2)
    })
    
    # Apply scaling if provided
    if (!is.null(scaling_vec)) {
      if (length(scaling_vec) == n_covs) {
        # Multiply only the median row (index 1) by the scaling vector
        for (k in 1:n_covs) {
          sum_vals[1, k, ] <- sum_vals[1, k, ] * scaling_vec[k]
        }
        subtitle <- paste0(subtitle, " (Scaled by SD)")
      } else {
        warning("Length of scaling_vec does not match number of covariates. Skipping scaling.")
      }
    }
    
    df_summary <- expand.grid(k = 1:n_covs, j = 1:n_cats)
    
    # Extract values safely
    df_summary$median <- apply(df_summary, 1, function(row) sum_vals[1, row["k"], row["j"]])
    df_summary$Credibility <- apply(df_summary, 1, function(row) sum_vals[2, row["k"], row["j"]])
    
    df_summary <- df_summary %>%
      mutate(
        Covariate = factor(cov_names[k], levels = rev(cov_names)),
        Category  = factor(cat_names[j], levels = cat_names)
      )
  } else if (is.data.frame(draws)) {
    # Handle summary table input (tab_coefs_5km style)
    df_summary <- draws
    
    # Try to map columns if standard names aren't present
    if (!"median" %in% names(df_summary) && "Estimate" %in% names(df_summary)) {
      df_summary$median <- df_summary$Estimate
    }
    
    if (!"Credibility" %in% names(df_summary)) {
      # Try to compute credibility from Upper/Lower or p-value if present
      if ("p_val" %in% names(df_summary)) {
        df_summary$Credibility <- abs(df_summary$p_val - 0.5) * 2
      } else {
        # Fallback to full size if we can't compute it
        df_summary$Credibility <- 1
      }
    }
    
    df_summary <- df_summary %>%
      mutate(
        Covariate = factor(Covariate, levels = rev(unique(Covariate))),
        Category  = factor(Category, levels = unique(Category))
      )
  } else {
    stop("Unsupported format for 'draws'. Expected a 3D array or a summary data frame.")
  }
  
  # Apply intercept filtering if requested
  if (exclude_intercept) {
    df_summary <- df_summary %>%
      filter(!grepl("intercept", Covariate, ignore.case = TRUE)) %>%
      mutate(
        Covariate = factor(Covariate, levels = rev(levels(Covariate)[levels(Covariate) %in% Covariate])),
        # Detect Neighborhood for faceting
        Panel = ifelse(grepl("Neighbor", Covariate), "Neighborhood Context", "Main Effects")
      )
  } else {
    df_summary <- df_summary %>%
      mutate(Panel = ifelse(grepl("Neighbor", Covariate), "Neighborhood Context", "Main Effects"))
  }
  
  # --- Apply custom sorting ---
  # 1. Sort Categories alphabetically, but force 'no_choice' to the end
  cats <- unique(as.character(df_summary$Category))
  cats_sorted <- sort(cats)
  if ("no_choice" %in% cats_sorted) {
    cats_sorted <- c(setdiff(cats_sorted, "no_choice"), "no_choice")
  }
  df_summary$Category <- factor(df_summary$Category, levels = cats_sorted)
  
  # 2. Sort Covariates: Main alphabetically, Neighborhood matching Categories
  covs <- unique(as.character(df_summary$Covariate))
  is_neighb <- grepl("Neighbor", covs, ignore.case = TRUE)
  
  neighb_covs <- covs[is_neighb]
  main_covs <- sort(covs[!is_neighb])
  
  # Target order for substring matching to ensure Y-axis matches X-axis
  no_choice_classes <- c("Waterbodies_marine", "Waterbodies_inland", "Wetlands_natural", "Natural_other", "NODATA")
  
  # Base the target order on cats_sorted so the neighborhood covariates exactly mirror the columns
  target_order <- c(cats_sorted, no_choice_classes)
  
  # Rank covariates based on which category they mention
  get_rank <- function(cov_name) {
    for (i in seq_along(target_order)) {
      if (grepl(target_order[i], cov_name, fixed = TRUE)) {
        return(i)
      }
    }
    return(999) # Fallback for unmatched
  }
  
  ranks <- sapply(neighb_covs, get_rank)
  neighb_covs_sorted <- neighb_covs[order(ranks, neighb_covs)]
  
  # ggplot y-axis draws from bottom to top, so to read top-to-bottom we reverse the levels
  cov_levels <- rev(c(main_covs, neighb_covs_sorted))
  df_summary$Covariate <- factor(df_summary$Covariate, levels = cov_levels)
  
  # Load ggnewscale for multiple legends
  require(ggnewscale)

  # Plot
  p <- ggplot(df_summary, aes(x = Category, y = Covariate)) +
    # Background tiles
    geom_tile(fill = "grey95", width = 1, height = 1) +
    
    # Main Effects tiles
    geom_tile(data = subset(df_summary, Panel == "Main Effects"),
              aes(width = Credibility, height = Credibility, fill = median)) +
    scale_fill_gradient2(low = "#7b3294", mid = "white", high = "#008837", 
                         midpoint = 0, name = "Est (Main)", 
                         guide = guide_colorbar(order = 1)) +
    
    # New scale for Neighborhood Context
    new_scale_fill() +
    geom_tile(data = subset(df_summary, Panel == "Neighborhood Context"),
              aes(width = Credibility, height = Credibility, fill = median)) +
    scale_fill_gradient2(low = "#7b3294", mid = "white", high = "#008837", 
                         midpoint = 0, name = "Est (Neighb.)", 
                         guide = guide_colorbar(order = 2)) +
                         
    # Borders
    geom_tile(fill = NA, color = "grey80", linewidth = 0.2) +
    
    # Facet with shared X but independent Y rows
    facet_grid(Panel ~ ., scales = "free_y", space = "free_y") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9, margin = margin(r = 2)),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      legend.position = "right",
      panel.spacing = unit(1, "lines")
    ) +
    labs(title = title,
         subtitle = subtitle,
         x = "Land-Use Category", 
         y = "Driver / Covariate")
  
  return(p)
}
