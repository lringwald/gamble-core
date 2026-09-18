library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(abind)

#' Advanced Random Effects Profiling for MNL
#' 
#' Creates a "Cloud" plot where each variable with random effects shows the distribution
#' of group-specific (e.g. Country) estimates against the EU-pooled mean.
#' 
#' @param res_combined A result list from mnlogit samplers containing postb and postb_pooled
#' @param re_idx Indices of covariates that have random effects
#' @param cov_names Character vector of covariate labels (understandable)
#' @param cat_names Character vector of land-use categories
#' @param group_names Character vector of group identifiers (e.g. NUTS0 codes)
#' @param title Plot title
#' @param sd_vec Optional vector of SDs for scaling if needed
#' 
#' @return A ggplot object
plot_mnl_hierarchical_profile <- function(res_combined, re_idx, cov_names, cat_names, group_names, 
                                          title = "Hierarchical Effect Profiles",
                                          sd_vec = NULL) {
  
  # 1. Extraction and Summary
  # res_combined$postb: [K x J-1 x G x Draws]
  # res_combined$postb_pooled: [K x J-1 x Draws]
  
  K <- dim(res_combined$postb)[1]
  J <- dim(res_combined$postb)[2]
  G <- dim(res_combined$postb)[3]
  
  # Validating labels
  if (is.null(cov_names)) cov_names <- dimnames(res_combined$postb)[[1]] %||% paste0("V", 1:K)
  if (is.null(cat_names)) cat_names <- dimnames(res_combined$postb)[[2]] %||% paste0("C", 1:J)
  if (is.null(group_names)) group_names <- as.character(1:G)
  
  # Build long-format summary of group effects
  eff_list <- list()
  for (k in re_idx) {
    for (j in 1:J) {
      # Group means
      g_means <- rowMeans(res_combined$postb[k, j, , ])
      
      # EU mean
      eu_mean <- mean(res_combined$postb_pooled[k, j, ])
      
      df_kj <- data.frame(
        Variable = cov_names[k],
        Category = cat_names[j],
        Group    = group_names,
        Mean     = g_means,
        EUMean   = eu_mean
      )
      eff_list[[length(eff_list) + 1]] <- df_kj
    }
  }
  
  df_long <- do.call(rbind, eff_list)
  
  # Assign Tiers for grouping/color
  tier_patterns <- c(
    Topography = "(Slope|Elevation|Aspect)",
    Economic   = "(GDP|Pop|Rent|Income|Price)",
    Policy     = "(allPA|GHM)",
    Neighborhood = "focal",
    Soil       = "(OC_TOP|ROO|AWC_TOP|VS)",
    Intercept  = "intercept"
  )
  
  df_long$Tier <- "Other"
  for (t in names(tier_patterns)) {
    df_long$Tier[grepl(tier_patterns[t], df_long$Variable, ignore.case = TRUE)] <- t
  }
  df_long$Tier <- factor(df_long$Tier, levels = c("Intercept", "Economic", "Topography", "Policy", "Neighborhood", "Soil", "Other"))
  
  # Ordering
  df_long$Variable <- factor(df_long$Variable, levels = rev(unique(df_long$Variable)))
  
  # Plot
  p <- ggplot(df_long) +
    # Reference EU-means as a background "shadow" point or bar
    geom_point(aes(x = EUMean, y = Variable), shape = 124, size = 4, color = "black", alpha = 0.5) +
    # Jittered group points
    geom_jitter(aes(x = Mean, y = Variable, color = Tier, group = Group), 
               height = 0.2, width = 0, alpha = 0.6, size = 1.5) +
    # Vertical zero line
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    facet_wrap(~Category, scales = "free_x") +
    scale_color_brewer(palette = "Set1") +
    theme_minimal(base_size = 11) +
    labs(title = title,
         subtitle = "Points = Country Estimates, Black tick = EU-Mean. Jittered vertically for visibility.",
         x = "Estimated Coefficient (Posterior Mean)", y = NULL, color = "Driver Group") +
    theme(panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 9),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
  
  return(p)
}

#' Random Effects Deviation Heatmap
#' 
#' Visualizes which groups (e.g. Countries) deviate most from the EU-average 
#' across all hierarchical variables.
#' 
#' @param res_combined A result list from mnlogit samplers
#' @param re_idx Indices of covariates with RE
#' @param cov_names Labels
#' @param cat_names Labels
#' @param group_names Labels
#' 
#' @return A ggplot object
plot_mnl_re_deviation_heatmap <- function(res_combined, re_idx, cov_names, cat_names, group_names, 
                                          title = "Response Deviations by Country") {
  
  K <- dim(res_combined$postb)[1]
  J <- dim(res_combined$postb)[2]
  G <- dim(res_combined$postb)[3]
  
  # Extract deviations
  dev_list <- list()
  for (k in re_idx) {
    for (j in 1:J) {
      g_means <- rowMeans(res_combined$postb[k, j, , ])
      eu_mean <- mean(res_combined$postb_pooled[k, j, ])
      
      df_kj <- data.frame(
        Variable  = cov_names[k],
        Category  = cat_names[j],
        Country   = group_names,
        Deviation = g_means - eu_mean
      )
      dev_list[[length(dev_list) + 1]] <- df_kj
    }
  }
  
  df_dev <- do.call(rbind, dev_list)
  
  # Reorder to keep categories and variables grouped
  df_dev$Variable <- factor(df_dev$Variable, levels = unique(df_dev$Variable))
  df_dev$Category <- factor(df_dev$Category, levels = unique(df_dev$Category))
  
  p <- ggplot(df_dev, aes(x = Variable, y = Country, fill = Deviation)) +
    geom_tile() +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "red4", 
                         midpoint = 0, name = "\u0394 Mean") +
    facet_wrap(~Category) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          axis.title = element_blank()) +
    labs(title = title,
         subtitle = "Red = Country response stronger/higher than EU-avg, Blue = Weaker/lower.")
  
  return(p)
}
