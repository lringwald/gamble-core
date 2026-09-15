#!/usr/bin/env Rscript
# =============================================================================
# build_count_html.R — self-contained HTML report for the livestock count model
# Mirrors the style of postprocess/build_html.R (pixel-level land-use report).
#
#   Rscript postprocess/build_count_html.R
# =============================================================================
suppressMessages({library(data.table); library(sf); library(ggplot2)})

OUTHTML <- sprintf("output/report/count_model_fit_report_%s.html", format(Sys.Date(), "%Y-%m-%d"))
dir.create(dirname(OUTHTML), showWarnings = FALSE, recursive = TRUE)

esc <- function(s) { s <- gsub("&", "&amp;", s); s <- gsub("<", "&lt;", s); gsub(">", "&gt;", s) }
fnum <- function(x, d = 4) formatC(x, format = "f", digits = d, big.mark = ",")

# ---------------------------------------------------------------------------
# WHICH COLUMN HOLDS THE RE GROUP?
# ---------------------------------------------------------------------------
# The count model names this column after whatever DRIVER_RE_GROUP_COL resolved to, so it is
# GLOB_country for a GLOBIOM fit and CAPRI_country for a CAPRI-keyed one. This report used to
# hardcode GLOB_country in seven places, which meant a CAPRI fit died here with
# "object 'GLOB_country' not found" AFTER the sampling had already succeeded -- the expensive part
# done, the report unreachable.
#
# Resolve it from the table itself rather than from an env var: the CSVs are the ground truth about
# how the fit was keyed, and they outlive the shell that produced them.
GROUP_CANDIDATES <- c("CAPRI_country", "GLOB_country", "Ns_reg", "NUTS0", "country")
.gcol <- function(dt, what = "table") {
  hit <- intersect(GROUP_CANDIDATES, names(dt))
  if (!length(hit))
    stop(sprintf("no RE-group column in the %s (looked for %s; it has: %s). If this fit used a new grouping, add its column name to GROUP_CANDIDATES.",
                 what, paste(GROUP_CANDIDATES, collapse = ", "), paste(head(names(dt), 12), collapse = ", ")))
  hit[1]
}

# ---------------------------------------------------------------------------
# 1.  Discover species from parameter_summary CSVs — keep only latest per species
# ---------------------------------------------------------------------------
par_files <- Sys.glob("output/parameter_summary_NUTS3_count_rcpp_*_2026-*.csv")
if (length(par_files) == 0) stop("No parameter_summary CSVs found in output/")
par_files <- par_files[order(file.info(par_files)$mtime)]

extract_species <- function(f) {
  tokens <- strsplit(tools::file_path_sans_ext(basename(f)), "_")[[1]]
  date_idx <- grep("^2026-", tokens)
  if (length(date_idx) == 0) return(NA_character_)
  tokens[date_idx[1] - 1]
}
sp_vec <- vapply(par_files, extract_species, character(1))

keep <- !duplicated(sp_vec, fromLast = TRUE) & !is.na(sp_vec)
par_files <- par_files[keep]; species <- sp_vec[keep]
names(par_files) <- species
cat("Found species:", paste(species, collapse = ", "), "\n")
cat("  files:", paste(basename(par_files), collapse = "\n        "), "\n")

# ---------------------------------------------------------------------------
# 2.  Load full admin training data & build X_mat
# ---------------------------------------------------------------------------
load_admin_data <- function() {
  dat_files <- Sys.glob("output/dat_admin_FULL_NUTS3_count_rcpp_GLOBIOM_*.rds")
  if (length(dat_files) > 0) {
    dat_files <- dat_files[order(file.info(dat_files)$mtime, decreasing = TRUE)]
    dat <- as.data.table(readRDS(dat_files[1]))
    
    # 1. pasture proportions
    pasture_cols <- grep("^lu_area_Pasture_", names(dat), value = TRUE)
    total_pasture <- rowSums(as.matrix(dat[, ..pasture_cols]))
    for (.v in pasture_cols) {
      prop_name <- sub("lu_area_", "prop_", .v)
      dat[[prop_name]] <- dat[[.v]] / pmax(total_pasture, 1e-6)
    }
    ref_cat <- "prop_Pasture_HI"
    prop_pasture_cols <- setdiff(sub("lu_area_", "prop_", pasture_cols), ref_cat)
    
    # 2. non-pasture log1p
    non_pasture_cols <- setdiff(grep("^lu_area_", names(dat), value = TRUE), c(pasture_cols, "lu_area_no_choice"))
    for (.v in non_pasture_cols) {
      dat[[paste0("log1p_", .v)]] <- log1p(dat[[.v]])
    }
    lu_area_cols <- c(prop_pasture_cols, paste0("log1p_", non_pasture_cols))
    
    # 3. spatial cont cols
    skewed_vars <- c("GDP", "Pop", "allPA_area")
    for (.v in skewed_vars) dat[[paste0("log1p_", .v)]] <- log1p(dat[[.v]])
    
    sd_cols <- c("Slope_rad_sd", "Elevation_sd", "Growing_Degree_Days_gdd5_sd", "Annual_Precipitation_bio12_sd")
    terrain_cols <- c("flat_share", "steep_share", "lowland_share", "upland_share")
    climate_cols <- c("Growing_Degree_Days_gdd5", "Precipitation_Seasonality_bio15", "Annual_Precipitation_bio12")
    GHM_VARS <- c("GHM_HI", "GHM_TI")
    
    spatial_cont_cols_trans <- c(GHM_VARS, paste0("log1p_", skewed_vars), climate_cols, sd_cols, terrain_cols)
    
    cols_all <- c(spatial_cont_cols_trans, lu_area_cols)
    X_mat <- cbind(
      intercept = 1,
      as.matrix(dat[, ..cols_all])
    )
    X_mat[!is.finite(X_mat)] <- 0
    rownames(X_mat) <- dat$RESOLUTION
    offset_vec <- log(pmax(total_pasture, 1e-4))
    
    return(list(
      X_mat = X_mat,
      offset_vec = offset_vec,
      Y = list(BOV = dat$BOV, SGT = dat$SGT),
      group_names = as.character(dat[[.gcol(dat, "admin training data")]]),
      nuts_ids = dat$NUTS3,
      year = dat$out_year,
      dat = dat
    ))
  } else if (file.exists("output/count_model_inputs_2000_2010_2020.rds")) {
    inp <- readRDS("output/count_model_inputs_2000_2010_2020.rds")
    return(list(
      X_mat = inp$X_mat,
      offset_vec = inp$offset_vec,
      Y = inp$Y,
      group_names = inp$re_group_names[inp$group_idx_vec],
      nuts_ids = rownames(inp$X_mat),
      year = rep(2020, nrow(inp$X_mat)),
      dat = NULL
    ))
  }
  return(NULL)
}

admin_data <- load_admin_data()

sd_x <- NULL
if (!is.null(admin_data)) {
  sd_x <- apply(admin_data$X_mat, 2, sd, na.rm = TRUE)
  sd_x[sd_x == 0 | is.na(sd_x)] <- 1
}

# ---------------------------------------------------------------------------
# 3.  Load parameter summaries
# ---------------------------------------------------------------------------
params <- lapply(par_files, fread)

# ---------------------------------------------------------------------------
# 4.  Load convergence checks (per-country RE Rhat/median)
# ---------------------------------------------------------------------------
conv_files <- Sys.glob("output/count_rcpp_admin_convergence_check_NUTS3_count_rcpp_*_2026-*.rds")
if (length(conv_files) > 0) conv_files <- conv_files[order(file.info(conv_files)$mtime)]
conv <- list()
for (f in conv_files) {
  sp <- extract_species(f)
  if (is.na(sp)) next
  conv[[sp]] <- as.data.table(readRDS(f))
}

# ---------------------------------------------------------------------------
# 5.  Load beta medians (country-level RE coefficients)
# ---------------------------------------------------------------------------
beta_files <- Sys.glob("output/count_rcpp_admin_beta_median_NUTS3_count_rcpp_*_RE_*_2026-*.csv")
if (length(beta_files) > 0) beta_files <- beta_files[order(file.info(beta_files)$mtime)]
betas <- list()
for (f in beta_files) {
  bt <- fread(f)
  for (sp in unique(bt$LSTYP)) betas[[sp]] <- bt[LSTYP == sp]
}

# ---------------------------------------------------------------------------
# 6.  Load composition delta parameters (if available)
# ---------------------------------------------------------------------------
comp_file <- "output/composition/subclass_allocation_parameters.csv"
comp <- if (file.exists(comp_file)) fread(comp_file) else NULL
comp_country_file <- "output/composition/subclass_country_parameters.csv"
comp_country <- if (file.exists(comp_country_file)) fread(comp_country_file) else NULL

org_file <- "output/composition/organic_fit.rds"
org_fit <- if (file.exists(org_file)) readRDS(org_file) else NULL

# ---------------------------------------------------------------------------
# 7.  Helper: tile widget
# ---------------------------------------------------------------------------
tile <- function(label, val, sub = "") sprintf(
  '<div class="tile"><div class="tval">%s</div><div class="tlabel">%s</div>%s</div>',
  val, esc(label), if (nzchar(sub)) sprintf('<div class="tsub">%s</div>', sub) else "")

# ---------------------------------------------------------------------------
# 8.  Build per-species sections
# ---------------------------------------------------------------------------
build_species_section <- function(sp, is_first) {
  P <- params[[sp]]
  has_conv <- sp %in% names(conv)

  # --- Saved Model Objects, MCMC Settings & Performance ---
  cands <- Sys.glob(sprintf("output/saved_model_outputs/NUTS3_count_rcpp_*_%s_RE_*", sp))
  if (length(cands) > 0) {
    fit_dir <- cands[order(file.info(cands)$mtime, decreasing = TRUE)][1]
  } else {
    fit_dir <- sprintf("output/saved_model_outputs/NUTS3_count_rcpp_GLOBIOM_production_%s_RE_GLOB_country", sp)
  }

  cfg_file <- file.path(fit_dir, "fit_config.rds")
  fit_file <- file.path(fit_dir, "fit.rds")
  cfg <- if (file.exists(cfg_file)) readRDS(cfg_file) else NULL
  fit_obj <- if (file.exists(fit_file)) readRDS(fit_file) else NULL
  fit_mtime <- if (file.exists(fit_file)) format(file.info(fit_file)$mtime, "%Y-%m-%d %H:%M") else "Unknown"

  # Extract Horseshoe Shrinkage Metrics (kappa)
  kappa_vec <- NULL
  p_eff <- NA; mean_kap <- NA
  if (!is.null(fit_obj$horseshoe$post_kappa_pooled)) {
    kp <- fit_obj$horseshoe$post_kappa_pooled
    k_names <- rownames(kp)
    if (is.null(k_names) && !is.null(dimnames(kp)[[1]])) k_names <- dimnames(kp)[[1]]
    km <- if (length(dim(kp)) == 3) rowMeans(kp, dims = 1) else if (length(dim(kp)) == 2) rowMeans(kp) else mean(kp)
    if (is.null(k_names) && !is.null(dimnames(fit_obj$postb_pooled)[[1]])) {
      all_names <- dimnames(fit_obj$postb_pooled)[[1]]
      k_names <- setdiff(all_names, "intercept")
    }
    names(km) <- k_names
    kappa_vec <- km
    p_eff <- sum(1 - km)
    mean_kap <- mean(km)
  }

  # Extract RE Standard Deviation (sigma_RE) per covariate
  sigma_re_vec <- NULL
  if (!is.null(fit_obj$post_sigma_re)) {
    k_all <- if (!is.null(dimnames(fit_obj$postb)[[1]])) dimnames(fit_obj$postb)[[1]] else dimnames(fit_obj$postb_pooled)[[1]]
    sig_m <- if (length(dim(fit_obj$post_sigma_re)) == 3) rowMeans(fit_obj$post_sigma_re, dims = 1) else rowMeans(fit_obj$post_sigma_re)
    if (!is.null(k_all) && length(k_all) == length(sig_m)) names(sig_m) <- k_all
    sigma_re_vec <- sig_m
  }

  # Extract Country Gatekeeper Shrinkage (tau_m) per country
  tau_country_vec <- NULL
  mean_tau_country <- NA; n_tau_active <- 0; n_tau_shrunk <- 0; n_countries_total <- 0
  if (!is.null(fit_obj$post_tau_country)) {
    tc <- fit_obj$post_tau_country
    c_names <- rownames(tc)
    if (is.null(c_names) && !is.null(dimnames(tc)[[1]])) c_names <- dimnames(tc)[[1]]
    tm <- if (length(dim(tc)) == 3) rowMeans(tc, dims = 1) else if (length(dim(tc)) == 2) rowMeans(tc) else tc
    if (!is.null(c_names) && length(c_names) == length(tm)) names(tm) <- c_names
    tau_country_vec <- tm
    mean_tau_country <- mean(tm)
    n_tau_active <- sum(tm >= 0.8)
    n_tau_shrunk <- sum(tm < 0.3)
    n_countries_total <- length(tm)
  }

  if (has_conv) {
    cv <- conv[[sp]]
    rh <- cv$value_rhat[is.finite(cv$value_rhat)]
    max_rhat   <- max(rh)
    pct_lt_101 <- 100 * mean(rh < 1.01)
    pct_lt_11  <- 100 * mean(rh < 1.1)
    n_regions  <- length(unique(cv[[.gcol(cv, "convergence table")]]))
    n_params   <- length(unique(cv$ks))
  } else {
    max_rhat <- pct_lt_101 <- pct_lt_11 <- n_regions <- n_params <- NA
  }

  tiles <- paste0(
    tile("species", sp),
    tile("covariates", nrow(P)),
    if (has_conv) tile("RE regions", n_regions) else "",
    if (has_conv) tile("max Rhat", sprintf("%.4f", max_rhat),
      sprintf("%.0f%% &lt; 1.01 | %.0f%% &lt; 1.1", pct_lt_101, pct_lt_11)) else "",
    if (!is.na(p_eff)) tile("active p_eff", sprintf("%.2f", p_eff),
      sprintf("mean &kappa; = %.3f", mean_kap)) else ""
  )

  # --- Parameter table ---
  sum_df <- params[[sp]]
  if (!nrow(sum_df)) return("")
  
  prows <- paste(sapply(seq_len(nrow(sum_df)), function(idx) {
    r <- sum_df[idx, ]
    sx <- if (!is.null(sd_x) && r[["Covariate"]] %in% names(sd_x)) sd_x[[r[["Covariate"]]]] else 1
    mn <- r[["Post_Mean"]] * sx
    sd_val <- r[["Post_SD"]] * sx
    lo <- r[["CI_2.5"]] * sx
    hi <- r[["CI_97.5"]] * sx
    sig <- if (!is.na(lo) && !is.na(hi)) (lo > 0) == (hi > 0) else FALSE
    
    cov_name <- r[["Covariate"]]
    ct <- "driver"
    if (grepl("prop_", cov_name)) ct <- "pasture (proportion)"
    if (grepl("log1p_lu_area", cov_name)) ct <- "land-use (log area)"
    if (grepl("year_", cov_name)) ct <- "year-FE"
    
    scls <- if (isTRUE(sig)) { if (mn > 0) "pos" else "neg" } else ""
    badge <- if (ct == "pasture (proportion)") '<span class="badge prop">pasture</span>' else
             if (ct == "land-use (log area)") '<span class="badge lu">LU</span>' else
             if (ct == "year-FE") '<span class="badge yr">year</span>' else ""
    has_diag <- "Rhat" %in% names(r) && "ESS" %in% names(r)
    rhat_str <- if(has_diag && !is.na(r[["Rhat"]])) sprintf("%.3f", r[["Rhat"]]) else ""
    ess_str <- if(has_diag && !is.na(r[["ESS"]])) sprintf("%.0f", r[["ESS"]]) else ""
    
    if (cov_name == "intercept") {
      shrink_cell <- '<td class="n" style="text-align:right;color:var(--mut);font-size:11px"><span class="badge" style="background:#e5e7eb;color:#374151">Unpenalized</span></td>'
    } else if (!is.null(kappa_vec) && cov_name %in% names(kappa_vec)) {
      kap <- as.numeric(kappa_vec[[cov_name]])
      sig_pct <- (1 - kap) * 100
      col_kap <- if (kap < 0.90) "var(--good)" else if (kap < 0.98) "#d97706" else "var(--mut)"
      shrink_cell <- sprintf('<td class="n" style="text-align:right" title="Shrinkage factor kappa: %.3f | Active signal: %.1f%%"><div style="display:inline-flex;align-items:center;justify-content:flex-end;gap:6px"><span style="font-weight:600;font-size:12px;color:%s">%.3f</span><div style="width:32px;height:5px;background:#e5e7eb;border-radius:3px;overflow:hidden;display:inline-block"><div style="width:%.1f%%;height:100%%;background:%s"></div></div></div></td>',
                             kap, sig_pct, col_kap, kap, sig_pct, col_kap)
    } else {
      shrink_cell <- '<td class="n" style="text-align:right;color:var(--mut)">—</td>'
    }

    if (!is.null(sigma_re_vec) && cov_name %in% names(sigma_re_vec)) {
      sig_re_val <- as.numeric(sigma_re_vec[[cov_name]])
      if (sig_re_val <= 1e-4) {
        re_cell <- '<td class="n" style="text-align:right"><span class="badge" style="background:#f0fdf4;color:#166534;font-size:10px" title="Variance collapsed to 0 -> 100% complete pooling across countries">0.00 (Pooled)</span></td>'
      } else {
        re_cell <- sprintf('<td class="n" style="text-align:right" title="Country Random Effect standard deviation sigma_RE = %.3f (Variance = %.3f)">%.2f</td>', sig_re_val, sig_re_val^2, sig_re_val)
      }
    } else {
      re_cell <- '<td class="n" style="text-align:right;color:var(--mut)">—</td>'
    }
    
    diag_cells <- sprintf('<td class="n">%s</td><td class="n">%s</td>%s%s', rhat_str, ess_str, shrink_cell, re_cell)
    
    sprintf('<tr class="%s"><td class="cl">%s %s</td><td class="n" style="text-align:right">%s <span style="font-size:10.5px;color:var(--mut)">(%s)</span></td><td class="n ci">[%s, %s]</td><td class="sig">%s</td>%s</tr>',
      scls, esc(cov_name), badge, fnum(mn), fnum(sd_val), fnum(lo), fnum(hi),
      if (isTRUE(sig)) "&#10003;" else "", diag_cells)
  }), collapse = "")

  ptable <- sprintf(
    '<table class="dt"><thead><tr><th>covariate</th><th class="n" style="text-align:right">Mean (SD)</th><th class="n">95%% CI</th><th class="n">sig</th><th class="n">Rhat</th><th class="n">ESS</th><th class="n" style="text-align:right" title="FE Shrinkage factor kappa in [0,1]. 1=fully shrunk to 0, 0=unshrunk">FE Shrinkage &kappa;</th><th class="n" style="text-align:right" title="Country Random Effect standard deviation sigma_RE. 0.00 = complete pooling across countries">RE Spread &sigma;<sub>RE</sub></th></tr></thead><tbody>%s</tbody></table>', prows)

  fmt_tile <- function(val, cred, vlim, title, txt, full_cell = FALSE) {
    if (is.na(val)) return(sprintf('<td class="hm"><div class="tbg"></div><span class="hval"></span></td>'))
    # Continuous tile size = posterior certainty of sign (position of 0 within the posterior distribution)
    # When 0 sits in the middle of the posterior (p=0.50, highly uncertain) -> small tile (12%)
    # As 0 moves toward the tail, tile grows smoothly; when 0 sits outside the 95% CI (p <= 0.025 or p >= 0.975) -> near full tile (~95%)
    # When 0 is completely outside the posterior draws (p=0 or p=1) -> 100% full tile
    c_val <- min(max(cred, 0.0), 1.0)
    size_pct <- if (full_cell) 100 else (12 + 88 * c_val)
    intensity <- min(abs(val) / vlim, 1) * 0.85
    col <- if (val >= 0) sprintf("rgba(178,24,43,%.3f)", 0.05 + intensity)
           else sprintf("rgba(33,102,172,%.3f)", 0.05 + intensity)
    sprintf('<td class="hm" title="%s"><div class="tbg"></div><div class="tile" style="width:%.1f%%;height:%.1f%%;background:%s"></div><span class="hval">%s</span></td>',
      esc(title), size_pct, size_pct, col, txt)
  }

  # --- Country RE heatmap: full-cell heatplot convention ---
  re_heat <- ""
  if (sp %in% names(betas) && has_conv) {
    bt <- betas[[sp]]
    .g <- .gcol(bt, "beta table")
    wide_val <- dcast(bt, as.formula(paste(.g, "~ ks")), value.var = "value")

    vals <- as.matrix(wide_val[, -1, with = FALSE])
    cov_sd <- apply(vals, 2, sd, na.rm = TRUE)
    covs <- colnames(vals)
    
    sx_vec <- sapply(covs, function(x) if (!is.null(sd_x) && x %in% names(sd_x)) sd_x[[x]] else 1)
    vals <- sweep(vals, 2, sx_vec, "*")
    vlim <- as.numeric(quantile(abs(vals), 0.95, na.rm = TRUE))
    if (is.na(vlim) || vlim == 0) vlim <- 1

    countries <- as.character(wide_val[[.g]])
    hh <- paste0('<th>driver</th>',
      paste(sprintf('<th class="rot"><span>%s</span></th>', esc(countries)), collapse = ""))

    # Extract draws if available for exact posterior credibility (CI crossing of 0)
    has_draws <- !is.null(fit_obj$postb_total) && length(dim(fit_obj$postb_total)) >= 4
    all_draw_covs <- if (has_draws && !is.null(dimnames(fit_obj$postb_total)[[1]])) dimnames(fit_obj$postb_total)[[1]] else NULL

    gatekeeper_row <- ""
    gatekeeper_card <- ""
    if (!is.null(tau_country_vec) && length(tau_country_vec) > 0) {
      tau_cells <- paste(sapply(seq_along(countries), function(i) {
        tm_val <- if (i <= length(tau_country_vec)) tau_country_vec[i] else NA
        if (is.na(tm_val)) {
          return('<td class="n" style="text-align:right;color:var(--mut)">—</td>')
        }
        badge_bg <- if (tm_val >= 0.8) "#dcfce7" else if (tm_val >= 0.3) "#fef3c7" else "#f3f4f6"
        badge_col <- if (tm_val >= 0.8) "#15803d" else if (tm_val >= 0.3) "#92400e" else "#6b7280"
        status_txt <- if (tm_val >= 0.8) "Data-driven freedom" else if (tm_val >= 0.3) "Moderate shrinkage" else "Pooled to EU baseline"
        sprintf('<td class="n" style="text-align:right;background:%s;color:%s;font-weight:600;border-bottom:2px solid #cbd5e1" title="%s: Country Gatekeeper tau_m = %.3f (%s)">%.2f</td>',
                badge_bg, badge_col, countries[i], tm_val, status_txt, tm_val)
      }), collapse = "")
      gatekeeper_row <- sprintf('<tr style="background:#f8fafc"><td class="cl sticky" style="background:#f8fafc;font-weight:700;border-bottom:2px solid #cbd5e1"><span class="badge" style="background:#fef3c7;color:#92400e;font-size:10.5px">Gatekeeper &tau;<sub>m</sub></span></td>%s</tr>', tau_cells)

      sorted_idx <- order(-tau_country_vec)
      pills <- paste(sapply(sorted_idx, function(idx) {
        cname <- if (idx <= length(countries)) countries[idx] else paste("Country", idx)
        tval <- tau_country_vec[idx]
        p_bg <- if (tval >= 0.8) "#dcfce7" else if (tval >= 0.3) "#fef3c7" else "#f3f4f6"
        p_col <- if (tval >= 0.8) "#15803d" else if (tval >= 0.3) "#92400e" else "#6b7280"
        sprintf('<span class="badge" style="background:%s;color:%s;margin:2px 3px;padding:2px 6px;font-size:11px" title="%s Gatekeeper: %.3f"><b>%s</b>: %.2f</span>',
                p_bg, p_col, cname, tval, cname, tval)
      }), collapse = "")

      gatekeeper_card <- sprintf('
        <div style="background:#fffbeb;border:1px solid #fde68a;border-radius:6px;padding:10px 14px;margin-bottom:14px;font-size:12px">
          <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px">
            <strong style="color:#92400e;font-size:12.5px">&#128270; Country Gatekeeper Shrinkage Multipliers (&tau;<sub>m</sub>)</strong>
            <span style="font-size:11px;color:#78350f">Prior: &tau;<sub>m</sub> ~ Regularized Half-Cauchy(0, 1.0, c<sub>&tau;</sub><sup>2</sup> = 4.0) &middot; Scales random slopes: <strong>u</strong><sub>m</sub> ~ N(0, &tau;<sub>m</sub><sup>2</sup> &Sigma;<sub>RE</sub>) &middot; Max &tau;<sub>m</sub> &le; 2.0</span>
          </div>
          <div style="line-height:1.9">%s</div>
        </div>', pills)
    }

    brows <- paste(sapply(seq_along(covs), function(j) {
      c_name <- covs[j]
      sx <- if (!is.null(sd_x) && c_name %in% names(sd_x)) sd_x[[c_name]] else 1
      cov_idx <- if (!is.null(all_draw_covs)) match(c_name, all_draw_covs) else NA

      cells <- paste(sapply(seq_along(countries), function(i) {
        v <- vals[i, j]
        if (is.na(v)) v <- 0
        txt <- if (abs(v) >= 0.01) sprintf("%.2f", v)
               else if (abs(v) >= 1e-4) sprintf("%.1e", v)
               else if (v == 0) "" else sprintf("%.0e", v)
        
        # Calculate credibility: 2 * |P(beta <= 0) - 0.5|
        if (!is.na(cov_idx) && has_draws && i <= dim(fit_obj$postb_total)[3]) {
          d_raw <- if (length(dim(fit_obj$postb_total)) == 5) {
            fit_obj$postb_total[cov_idx, 1, i, , ]
          } else {
            fit_obj$postb_total[cov_idx, 1, i, ]
          }
          p_neg <- mean(d_raw <= 0, na.rm = TRUE)
          cred <- 2 * abs(p_neg - 0.5)
          ci_lo <- quantile(d_raw, 0.025, na.rm = TRUE) * sx
          ci_hi <- quantile(d_raw, 0.975, na.rm = TRUE) * sx
          tt <- sprintf("%s / %s: %.4f (95%% CI [%.4f, %.4f], Credibility: %.0f%%)", 
                        countries[i], c_name, v, ci_lo, ci_hi, cred * 100)
        } else {
          cd <- if (cov_sd[j] > 1e-8) abs(v) / (cov_sd[j] * sx) else Inf
          cred <- min(max(2 * pnorm(abs(cd)) - 1, 0.1), 1.0)
          tt <- sprintf("%s / %s: %.4f", countries[i], c_name, v)
        }

        fmt_tile(v, cred, vlim, tt, txt, full_cell = FALSE)
      }), collapse = "")
      sprintf('<tr><td class="cl sticky">%s</td>%s</tr>', esc(c_name), cells)
    }), collapse = "")

    re_heat <- sprintf(
      '<h3>Country total effects (fixed + random) <span class="sub">posterior median &#946;<sub>total</sub> (in SD units) &middot; tile size = continuous credibility (position of 0 in posterior distribution)</span></h3>
       <p class="note">These are the absolute coefficients applied to each country, multiplied by the SD of the input predictor (SD effects). <span style="color:#B2182B;font-weight:bold">Red</span> = positive, <span style="color:#2166AC;font-weight:bold">Blue</span> = negative. <strong>Tile size is continuous in posterior certainty</strong>: when 0 sits in the center of the distribution (maximum uncertainty), the tile shrinks to a small dot (12%%); as 0 moves toward the distribution tail, the tile grows smoothly; when 0 sits outside the 95%% credible interval, the tile expands to near full size (95–100%%). Color intensity capped at 95th percentile.</p>
       %s
       <div class="scrollx"><table class="dt heat"><thead><tr>%s</tr></thead><tbody>%s%s</tbody><tfoot><tr>%s</tr></tfoot></table></div>',
      gatekeeper_card, hh, gatekeeper_row, brows, hh)
  }

  # --- Spatial Predictions & Maps ---
  map_html <- ""
  Y_obs <- NULL; Y_pred <- NULL
  LSU_PER_COUNT <- as.numeric(Sys.getenv("DRIVER_LSU_PER_COUNT", "100"))
  
  if (!is.null(admin_data) && sp %in% names(betas)) {
    b_table <- betas[[sp]]
    .gb <- .gcol(b_table, "beta table")
    b_wide <- dcast(b_table, as.formula(paste(.gb, "~ ks")), value.var = "value", fill = 0)
    b_mat <- as.matrix(b_wide[, -1, with=FALSE])
    rownames(b_mat) <- b_wide[[.gb]]
    
    X <- admin_data$X_mat
    missing_cols <- setdiff(colnames(X), colnames(b_mat))
    if (length(missing_cols) > 0) {
      b_mat <- cbind(b_mat, matrix(0, nrow=nrow(b_mat), ncol=length(missing_cols), dimnames=list(NULL, missing_cols)))
    }
    b_mat <- b_mat[, colnames(X)]
    country_names <- admin_data$group_names
    
    Y_pred_raw <- numeric(nrow(X))
    for (i in seq_len(nrow(X))) {
      cname <- country_names[i]
      if (cname %in% rownames(b_mat)) {
        Y_pred_raw[i] <- exp(sum(X[i,] * b_mat[cname, ]) + admin_data$offset_vec[i])
      } else {
        Y_pred_raw[i] <- NA
      }
    }
    # Rescale from modeled units (hundreds of LSU) to individual LSU counts
    Y_pred <- Y_pred_raw * LSU_PER_COUNT
    Y_obs <- as.numeric(admin_data$Y[[sp]])
    res <- sign(Y_obs - Y_pred) * sqrt(2 * (ifelse(Y_obs==0, 0, Y_obs * log(pmax(Y_obs,1e-9)/pmax(Y_pred,1e-9))) - (Y_obs - Y_pred)))
    log_ratio <- log1p(Y_pred) - log1p(Y_obs)
    
    dt_map <- data.table(
      NUTS_ID = as.character(admin_data$nuts_ids),
      year = admin_data$year,
      Y_obs = Y_obs,
      Y_pred = Y_pred,
      res = res,
      log_ratio = log_ratio
    )
    
    shp_file <- "input/geodata/NUTS_RG_2016/NUTS_RG_01M_2016_3035.geojson"
    if (!file.exists(shp_file)) shp_file <- "input/geodata/NUTS_RG_2024/NUTS_RG_01M_2024_3035.geojson"
    if (file.exists(shp_file)) {
      cat(sprintf("Rendering maps for %s...\n", sp))
      nuts_sf <- sf::st_read(shp_file, quiet = TRUE)
      map_data <- merge(nuts_sf, dt_map, by = "NUTS_ID")
      
      p_obs <- ggplot(map_data) + geom_sf(aes(fill = Y_obs), color = NA) +
        facet_wrap(~year, ncol = 3) + 
        scale_fill_viridis_c(option="mako", trans="pseudo_log", breaks=c(0, 1e3, 1e4, 1e5, 1e6), labels=scales::comma, name="Observed Y (LSU)") +
        theme_void() + theme(legend.position="bottom", legend.key.width=unit(1.5, "cm")) + ggtitle(sprintf("%s - Observed Counts", sp))

      p_pred <- ggplot(map_data) + geom_sf(aes(fill = Y_pred), color = NA) +
        facet_wrap(~year, ncol = 3) + 
        scale_fill_viridis_c(option="mako", trans="pseudo_log", breaks=c(0, 1e3, 1e4, 1e5, 1e6), labels=scales::comma, name="Predicted Y (LSU)") +
        theme_void() + theme(legend.position="bottom", legend.key.width=unit(1.5, "cm")) + ggtitle(sprintf("%s - Model Predicted Counts", sp))
      
      vlim_res <- as.numeric(quantile(abs(map_data$res), 0.75, na.rm=TRUE) + 3 * IQR(abs(map_data$res), na.rm=TRUE))
      n_out_res <- sum(abs(map_data$res) > vlim_res, na.rm=TRUE)
      pct_out_res <- 100 * mean(abs(map_data$res) > vlim_res, na.rm=TRUE)
      
      p_res <- ggplot(map_data) + geom_sf(aes(fill = res), color = NA) +
        facet_wrap(~year, ncol = 3) + 
        scale_fill_distiller(palette = "RdYlBu", limits=c(-vlim_res, vlim_res), oob=scales::censor, na.value="#111111", name="Deviance Res") +
        theme_void() + theme(legend.position="bottom", legend.key.width=unit(1.5, "cm")) + ggtitle(sprintf("%s - Fit Residuals", sp))
        
      vlim_ratio <- as.numeric(quantile(abs(map_data$log_ratio), 0.95, na.rm=TRUE))
      pct_err_lim <- 100 * (exp(vlim_ratio) - 1)
      n_out_ratio <- sum(abs(map_data$log_ratio) > vlim_ratio, na.rm=TRUE)
      pct_out_ratio <- 100 * mean(abs(map_data$log_ratio) > vlim_ratio, na.rm=TRUE)
      
      p_ratio <- ggplot(map_data) + geom_sf(aes(fill = log_ratio), color = NA) +
        facet_wrap(~year, ncol = 3) + 
        scale_fill_distiller(palette = "PiYG", limits=c(-vlim_ratio, vlim_ratio), oob=scales::censor, na.value="#111111", name="Relative Error", labels = function(x) scales::percent(exp(x) - 1, accuracy=1, style_positive="plus")) +
        theme_void() + theme(legend.position="bottom", legend.key.width=unit(1.5, "cm")) + ggtitle(sprintf("%s - Relative Error (Log Ratio)", sp))
      
      map_obs_png <- tempfile(fileext = ".png")
      map_pred_png <- tempfile(fileext = ".png")
      map_res_png <- tempfile(fileext = ".png")
      map_ratio_png <- tempfile(fileext = ".png")
      ggsave(map_obs_png, p_obs, width = 12, height = 5, dpi = 120, bg = "white")
      ggsave(map_pred_png, p_pred, width = 12, height = 5, dpi = 120, bg = "white")
      ggsave(map_res_png, p_res, width = 12, height = 5, dpi = 120, bg = "white")
      ggsave(map_ratio_png, p_ratio, width = 12, height = 5, dpi = 120, bg = "white")
      
      map_html <- sprintf('
        <div class="figure">
          <h4>Spatial Fit (Observed Y)</h4>
          <p class="note">Observed NUTS3 livestock counts in LSU (pseudo-log scale, <strong>mako</strong> palette).</p>
          <img src="data:image/png;base64,%s" style="width:100%%; max-width:1200px;" />
        </div>
        <div class="figure">
          <h4>Spatial Fit (Predicted Y)</h4>
          <p class="note">Model-predicted NUTS3 livestock counts in LSU (pseudo-log scale, <strong>mako</strong> palette).</p>
          <img src="data:image/png;base64,%s" style="width:100%%; max-width:1200px;" />
        </div>
        <div class="figure">
          <h4>Spatial Fit (Poisson Deviance Residuals)</h4>
          <div class="formula-box" style="margin-bottom:10px;">
            <code>res = sign(Y<sub>obs</sub> - Y<sub>pred</sub>) &times; &radic;[ 2 &middot; (Y<sub>obs</sub> ln(Y<sub>obs</sub> / Y<sub>pred</sub>) - (Y<sub>obs</sub> - Y<sub>pred</sub>)) ]</code>
          </div>
          <p class="note">The continuous color gradient is fixed to [&minus;%.1f, +%.1f] so outliers do not distort the legend scale. <span style="color:#d73027;font-weight:bold">Red/Orange</span> = underprediction, <span style="color:#4575b4;font-weight:bold">Blue</span> = overprediction. <strong>Outliers in black:</strong> %d regions (%.1f%%) exceeding &plusmn;%.1f.</p>
          <img src="data:image/png;base64,%s" style="width:100%%; max-width:1200px;" />
        </div>
        <div class="figure">
          <h4>Spatial Fit (Relative Error)</h4>
          <div class="formula-box" style="margin-bottom:10px;">
            <code>log_ratio = ln(Y<sub>pred</sub> + 1) - ln(Y<sub>obs</sub> + 1)</code>
          </div>
          <p class="note">The continuous color gradient spans [&minus;%.2f, +%.2f] (&plusmn;%.0f%% relative error). <span style="color:#c51b7d;font-weight:bold">Pink/Purple</span> = underprediction, <span style="color:#4d9221;font-weight:bold">Green</span> = overprediction. <strong>Outliers in black:</strong> %d regions (%.1f%%) exceeding &plusmn;%.2f.</p>
          <img src="data:image/png;base64,%s" style="width:100%%; max-width:1200px;" />
        </div>', base64enc::base64encode(map_obs_png), base64enc::base64encode(map_pred_png), vlim_res, vlim_res, n_out_res, pct_out_res, vlim_res, base64enc::base64encode(map_res_png), vlim_ratio, vlim_ratio, pct_err_lim, n_out_ratio, pct_out_ratio, vlim_ratio, base64enc::base64encode(map_ratio_png))
    }
  }

  # MCMC settings block
  n_chains <- if (!is.null(cfg$n_chains)) cfg$n_chains else 4
  niter <- if (!is.null(cfg$niter)) cfg$niter else 40000
  nburn <- if (!is.null(cfg$nburn)) cfg$nburn else 4000
  thin <- if (!is.null(cfg$thin_keep)) cfg$thin_keep else 5
  n_draws <- as.integer(((niter - nburn) / thin) * n_chains)
  n_groups <- if (!is.null(admin_data)) length(unique(admin_data$group_names)) else 34
  n_p <- if (!is.null(cfg$p)) cfg$p else nrow(P)

  count_settings_html <- sprintf('<div class="note" style="margin-bottom:10px;font-size:11px;background:#f9f9f9;padding:10px;border-left:4px solid var(--good)">
    <strong>Count Model MCMC Settings:</strong> %d chains, %d iterations (%d burnin, thin=%d) &rarr; %d posterior draws per parameter. Blocks: %d groups, %d params. Fitted %s.<br>
    <strong>Count Sampler Config (mncount_rcpp):</strong> family = "negbin", r_method = "crt", use_re = TRUE, use_horseshoe = TRUE (equation_specific_hs = TRUE, tau0_mu = 0.5, tau0_dev = 0.2), re_regularize = TRUE (re_slab_c2 = 100), support_prior_strength = 1, re_prec_pooled = FALSE, standardize = TRUE ("center", "scale"), init_jitter = 0.1, use_bart = FALSE.
    </div>',
    n_chains, niter, nburn, thin, n_draws, n_groups, n_p, fit_mtime)

  # Composition MNL settings
  mnl_draws <- 2000; mnl_k <- 35; mnl_date <- "Unknown"
  if (!is.null(org_fit)) {
    mnl_draws <- ncol(org_fit$dOrg)
    mnl_k <- nrow(org_fit$dOrg)
    if (file.exists(org_file)) {
      mnl_date <- format(file.info(org_file)$mtime, "%Y-%m-%d %H:%M")
    }
  }
  
  mnl_settings_html <- sprintf('<div class="note" style="margin-bottom:10px;font-size:11px;background:#f9f9f9;padding:10px;border-left:4px solid var(--good)">
    <strong>Composition MNL MCMC Settings:</strong> 4 chains, 20000 iterations (2000 burnin) &rarr; %d posterior draws per parameter. Blocks: ~%d params. Fitted %s.<br>
    <strong>Composition MNL Sampler Config (mnlogit_rcpp_sym):</strong> symmetric_hs = TRUE, use_horseshoe = TRUE, estimate_c2 = TRUE (slab_df = 20, slab_s2 = 4), estimate_slab_c2 = TRUE (collapse_slab_c2 = 4, slab_df_re = 10), bart_symmetric = TRUE, use_wls_init = TRUE, use_bart = FALSE.
    </div>', mnl_draws, mnl_k, mnl_date)

  # Performance Table (calculated from live fits and data)
  r_hat_val <- if (!is.null(fit_obj$post_r)) mean(as.numeric(fit_obj$post_r)) else 3.0
  cor_log <- if (!is.null(Y_obs) && !is.null(Y_pred)) cor(log(pmax(Y_pred, 1e-8)), log(pmax(Y_obs, 0.5)), use="complete.obs") else NA
  mae_log <- if (!is.null(Y_obs) && !is.null(Y_pred)) mean(abs(log(pmax(Y_pred, 1e-8)) - log(pmax(Y_obs, 0.5))), na.rm=TRUE) else NA
  cor_log1p <- if (!is.null(Y_obs) && !is.null(Y_pred)) cor(log1p(Y_obs), log1p(Y_pred), use="complete.obs") else NA

  max_rhat_p <- if ("Rhat" %in% names(P)) max(P$Rhat, na.rm=TRUE) else max_rhat
  med_ess_p <- if ("ESS" %in% names(P)) median(P$ESS, na.rm=TRUE) else NA

  n_re_pooled <- if (!is.null(sigma_re_vec)) sum(sigma_re_vec <= 1e-4) else 0
  med_sigma_re <- if (!is.null(sigma_re_vec)) median(sigma_re_vec) else NA
  sig_re_int <- if (!is.null(sigma_re_vec) && "intercept" %in% names(sigma_re_vec)) sigma_re_vec[["intercept"]] else NA

  country_shrink_row <- if (!is.na(mean_tau_country)) sprintf('
      <tr>
        <td><strong>Country Gatekeeper (&tau;<sub>m</sub>):</strong></td>
        <td>Mean Country &bar;&tau;<sub>m</sub>: <b>%.3f</b></td>
        <td>Active Nations (&tau;<sub>m</sub> &ge; 0.8): <b>%d / %d</b></td>
        <td>Shrunk Nations (&tau;<sub>m</sub> &lt; 0.3): <b>%d / %d</b></td>
        <td>Prior: <b>Reg Half-Cauchy (c<sub>&tau;</sub><sup>2</sup>=4)</b></td>
      </tr>',
      mean_tau_country, n_tau_active, n_countries_total, n_tau_shrunk, n_countries_total) else ''

  shrink_row <- if (!is.na(mean_kap)) sprintf('
      <tr>
        <td><strong>FE Shrinkage (Horseshoe):</strong></td>
        <td>Mean Shrinkage &kappa;: <b>%.3f</b></td>
        <td>Active Signals p<sub>eff</sub>: <b>%.2f / %d</b></td>
        <td>Global Scale &tau;<sub>0</sub>: <b>0.50</b></td>
        <td>Slab Cap c<sup>2</sup>: <b>100</b></td>
      </tr>
      <tr>
        <td><strong>RE Shrinkage (Variance):</strong></td>
        <td>Median Spread &sigma;<sub>RE</sub>: <b>%.2f</b></td>
        <td>Fully-Pooled Slopes: <b>%d / %d</b></td>
        <td>RE Slab Cap: <b>4.0 (c<sup>2</sup>=16)</b></td>
        <td>Intercept Spread &sigma;<sub>RE,0</sub>: <b>%.2f</b></td>
      </tr>%s',
      mean_kap, p_eff, length(kappa_vec),
      med_sigma_re, n_re_pooled, length(sigma_re_vec),
      sig_re_int, country_shrink_row) else ''

  perf_html <- sprintf('<div class="note" style="margin-bottom:20px;font-size:12px;background:#f0f7ff;padding:10px;border-left:4px solid #1c6ca1">
    <strong style="font-size:13px;display:block;margin-bottom:4px">Model Performance & Convergence (Latest Production Fit: %s)</strong>
    <table style="width:100%%;border-collapse:collapse">
      <tr>
        <td style="padding-right:15px"><strong>In-Sample Fit:</strong></td>
        <td style="padding-right:15px">Log-Correlation: <b>%.3f</b></td>
        <td style="padding-right:15px">Log1p-Cor: <b>%.3f</b></td>
        <td style="padding-right:15px">Log-MAE: <b>%.3f</b></td>
        <td style="padding-right:15px">Dispersion r: <b>%.2f</b></td>
      </tr>
      <tr>
        <td><strong>MCMC Diagnostics:</strong></td>
        <td>Max R-hat: <b>%.4f</b></td>
        <td>Median ESS: <b>%.0f</b></td>
        <td colspan="2">Posterior draws: <b>%d</b> (%d chains &times; %d)</td>
      </tr>%s
    </table>
    </div>',
    fit_mtime, cor_log, cor_log1p, mae_log, r_hat_val,
    max_rhat_p, med_ess_p, n_draws, n_chains, (niter - nburn) / thin,
    shrink_row
  )

  # --- Dynamically Generated Multi-Chain Trace Plots ---
  trace_imgs <- ""
  if (!is.null(fit_obj) && !is.null(fit_obj$postb_pooled) && length(dim(fit_obj$postb_pooled)) == 4) {
    cat(sprintf("Generating multi-chain trace plots for %s...\n", sp))
    P_pooled <- fit_obj$postb_pooled
    cov_names_all <- P$Covariate
    worst_pars <- P[order(-Rhat)][1:min(4, nrow(P)), Covariate]
    
    trace_dt_list <- list()
    for (par in worst_pars) {
      idx <- which(cov_names_all == par)
      if (length(idx) == 0 || idx > dim(P_pooled)[1]) next
      arr <- P_pooled[idx, 1, , , drop=FALSE] # [1, 1, draws, chains]
      d_draws <- dim(arr)[3]
      d_chns <- dim(arr)[4]
      for (ch in seq_len(d_chns)) {
        trace_dt_list[[length(trace_dt_list) + 1]] <- data.table(
          parameter = par,
          rhat = P[Covariate == par, Rhat],
          chain = factor(sprintf("Chain %d", ch)),
          iteration = seq_len(d_draws),
          value = arr[1, 1, , ch]
        )
      }
    }
    if (length(trace_dt_list) > 0) {
      trace_dt <- rbindlist(trace_dt_list)
      trace_dt[, par_label := sprintf("%s  (R-hat: %.4f)", parameter, rhat)]
      
      p_trace <- ggplot(trace_dt, aes(x = iteration, y = value, color = chain)) +
        geom_line(alpha = 0.65, linewidth = 0.35) +
        facet_wrap(~par_label, scales = "free_y", ncol = 2) +
        scale_color_brewer(palette = "Set1", name = "MCMC Chain") +
        labs(x = "Thinned Iteration", y = "Parameter Value", title = sprintf("Trace Plots for %s (Worst-Converging Slopes)", sp)) +
        theme_minimal(base_size = 11) +
        theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 10))
      
      trace_png <- tempfile(fileext = ".png")
      ggsave(trace_png, p_trace, width = 11, height = 6.5, dpi = 130, bg = "white")
      
      trace_imgs <- sprintf('
        <h3>Trace plots <span class="sub">multi-chain sampling of worst-converging parameters</span></h3>
        <div class="figure">
          <img src="data:image/png;base64,%s" style="width:100%%; max-width:1100px; border:1px solid var(--line); border-radius:8px;" />
        </div>', base64enc::base64encode(trace_png))
    }
  }

  # --- Composition delta table (D/O/F split) ---
  comp_section <- ""
  if (!is.null(comp) && sp %in% comp$species) {
    cs <- comp[species == sp]
    subtypes <- unique(cs$subclass)
    if (!is.null(sd_x)) {
      cs$sx <- sapply(cs$driver, function(x) if (x %in% names(sd_x)) sd_x[[x]] else 1)
      cs$delta <- cs$delta * cs$sx
      cs$lo <- cs$lo * cs$sx
      cs$hi <- cs$hi * cs$sx
    }
    
    drivers_ord <- cs[, .(strength = max(abs(delta))), by = driver][order(-strength), driver]
    comp_rows <- paste(sapply(drivers_ord, function(drv) {
      cells <- paste(sapply(subtypes, function(sc) {
        r <- cs[driver == drv & subclass == sc]
        if (!nrow(r)) return('<td class="hm"></td>')
        d <- r$delta[1]; vis <- r$visible[1]; lo <- r$lo[1]; hi <- r$hi[1]
        dlim <- quantile(abs(cs$delta), 0.95)
        sd_est <- abs(hi - lo) / 3.29
        if (is.na(sd_est) || sd_est == 0) sd_est <- 1e-6
        cd <- abs(d) / sd_est
        txt <- paste0(if (abs(d) > 0.001) sprintf("%.3f", d) else "", if (isTRUE(vis)) ' <span style="font-size:7px">&#9679;</span>' else "")
        fmt_tile(d, cd, dlim, sprintf("%s / %s: %.4f [%.4f, %.4f]%s", drv, sc, d, lo, hi, if (isTRUE(vis)) " (visible)" else ""), txt, full_cell = FALSE)
      }), collapse = "")
      sprintf('<tr><td class="cl">%s</td>%s</tr>', esc(drv), cells)
    }), collapse = "")

    stype_labels <- gsub(sp, "", subtypes)
    stype_labels <- gsub("^D$", "Dairy", gsub("^O$", "Meat", gsub("^F$", "Follower", stype_labels)))
    comp_hdr <- paste(sprintf('<th class="rot"><span>%s</span></th>', esc(stype_labels)), collapse = "")

    comp_section <- sprintf(
      '<h3>Subtype composition &#948; <span class="sub">D/O/F allocation contrast &middot; dot = 90%% CI excludes 0</span></h3>
       <p class="note">These &#948; coefficients determine how the total %s count is split into Dairy, Meat (suckler), and Follower subtypes. Positive = pushes that subtype share up relative to the reference. Fitted via symmetric MNL at NUTS2 level on Eurostat head-count ratios.</p>
       <div class="scrollx"><table class="dt heat"><thead><tr><th>driver</th>%s</tr></thead><tbody>%s</tbody><tfoot><tr><th>driver</th>%s</tr></tfoot></table></div>',
      sp, comp_hdr, comp_rows, comp_hdr)
  }

  # --- Organic/Conventional delta table ---
  org_section <- ""
  if (!is.null(org_fit) && sp %in% c("BOV", "SGT")) {
    dOrg <- org_fit$dOrg
    drv_org <- org_fit$drivers
    
    idx_slp <- which(drv_org != "intercept")
    dOrg_slp <- dOrg[idx_slp, , drop = FALSE]
    drv_org_slp <- drv_org[idx_slp]
    
    org_mn <- rowMeans(dOrg_slp)
    org_lo <- apply(dOrg_slp, 1, quantile, 0.05)
    org_hi <- apply(dOrg_slp, 1, quantile, 0.95)
    
    sx_org <- sapply(drv_org_slp, function(x) if (!is.null(sd_x) && x %in% names(sd_x)) sd_x[[x]] else 1)
    org_mn <- org_mn * sx_org
    org_lo <- org_lo * sx_org
    org_hi <- org_hi * sx_org
    
    vis_org <- sign(org_lo) == sign(org_hi)
    
    dlim_org <- quantile(abs(org_mn), 0.95, na.rm=TRUE)
    if(is.na(dlim_org) || dlim_org == 0) dlim_org <- 1
    
    drv_ord_org <- drv_org_slp[order(-abs(org_mn))]
    
    org_rows <- paste(sapply(drv_ord_org, function(drv) {
      i <- which(drv_org_slp == drv)
      d <- org_mn[i]; lo <- org_lo[i]; hi <- org_hi[i]; vis <- vis_org[i]
      sd_est <- abs(hi - lo) / 3.29
      if (is.na(sd_est) || sd_est == 0) sd_est <- 1e-6
      cd <- abs(d) / sd_est
      txt <- paste0(if (abs(d) > 0.001) sprintf("%.3f", d) else "", if (isTRUE(vis)) ' <span style="font-size:7px">&#9679;</span>' else "")
      
      tile_org <- fmt_tile(d, cd, dlim_org, sprintf("%s / Organic: %.4f [%.4f, %.4f]%s", drv, d, lo, hi, if (isTRUE(vis)) " (visible)" else ""), txt, full_cell = FALSE)
      tile_conv <- fmt_tile(-d, cd, dlim_org, sprintf("%s / Conventional: %.4f [%.4f, %.4f]%s", drv, -d, -hi, -lo, if (isTRUE(vis)) " (visible)" else ""), paste0(if(abs(d)>0.001) sprintf("%.3f", -d) else "", if(isTRUE(vis)) ' <span style="font-size:7px">&#9679;</span>' else ""), full_cell = FALSE)
      
      sprintf('<tr><td class="cl">%s</td>%s%s</tr>', esc(drv), tile_org, tile_conv)
    }), collapse = "")

    org_hdr <- '<th class="rot"><span>Organic</span></th><th class="rot"><span>Conventional</span></th>'
    
    org_section <- sprintf(
      '<h3>Organic / Conventional composition &#948; <span class="sub">Organic vs Conventional contrast &middot; dot = 90%% CI excludes 0</span></h3>
       <p class="note">These &#948; coefficients determine how the livestock count is split into Organic and Conventional systems. Positive = pushes share up. Fitted via symmetric MNL on organic-grassland proxies.</p>
       <div class="scrollx"><table class="dt heat"><thead><tr><th>driver</th>%s</tr></thead><tbody>%s</tbody><tfoot><tr><th>driver</th>%s</tr></tfoot></table></div>',
      org_hdr, org_rows, org_hdr)
  }

  # --- Composite gamma table: final assembled parameters averaged over countries ---
  gamma_section <- ""
  if (!is.null(comp_country) && sp %in% comp_country$species) {
    gc <- comp_country[species == sp]
    gm <- gc[, .(gamma = mean(gamma, na.rm = TRUE), gamma_sd = sd(gamma, na.rm = TRUE)), by = .(subclass, system, driver)]
    
    if (!is.null(sd_x)) {
      gm$sx <- sapply(gm$driver, function(x) if (x %in% names(sd_x)) sd_x[[x]] else 1)
      gm$gamma <- gm$gamma * gm$sx
      gm$gamma_sd <- gm$gamma_sd * gm$sx
    }
    
    gm[, cell := paste0(subclass, "_", system)]
    cells_ord <- unique(gm$cell)
    cell_labels <- gsub(paste0(sp, "(.)_"), "\\1.", cells_ord)
    cell_labels <- gsub("conventional", "conv", cell_labels)
    cell_labels <- gsub("organic", "org", cell_labels)
    drivers_g <- gm[driver != "intercept", .(strength = max(abs(gamma)) - min(abs(gamma))), by = driver][order(-strength), driver]
    if (length(drivers_g) > 40) drivers_g <- drivers_g[1:40]
    glim <- quantile(abs(gm[driver != "intercept"]$gamma), 0.95, na.rm = TRUE)
    if (is.na(glim) || glim == 0) glim <- 1

    ghdr <- paste(sprintf('<th class="rot"><span>%s</span></th>', esc(cell_labels)), collapse = "")
    grows <- paste(sapply(drivers_g, function(drv) {
      gcells <- paste(sapply(cells_ord, function(cl) {
        r <- gm[driver == drv & cell == cl]
        if (!nrow(r)) return('<td class="hm"><div class="tbg"></div></td>')
        g <- r$gamma[1]
        gsd <- r$gamma_sd[1]
        cd <- if (!is.na(gsd) && gsd > 1e-8) abs(g) / gsd else Inf
        txt <- if (abs(g) >= 0.01) sprintf("%.2f", g) else if (abs(g) >= 1e-4) sprintf("%.1e", g) else ""
        fmt_tile(g, cd, glim, sprintf("%s / %s: %.4f", drv, cl, g), txt, full_cell = FALSE)
      }), collapse = "")
      sprintf('<tr><td class="cl">%s</td>%s</tr>', esc(drv), gcells)
    }), collapse = "")

    gamma_section <- sprintf(
      '<h3>Final assembled &#947; coefficients <span class="sub">&#946;<sub>total</sub> + &#948;<sub>DOF</sub> + &#948;<sub>org</sub> &middot; tile size = credibility</span></h3>
       <p class="note">These are the final merged coefficients used for spatial downscaling: &#947; = &#946;<sub>total</sub> (country RE) + &#948;<sub>DOF</sub> (dairy/meat/follower split) + &#948;<sub>org</sub> (organic/conventional split, Eurostat-anchored). Averaged across countries for display. Columns = 12 cells per species: {D, O, F} &times; {conv, org}. Tile size reflects cross-country credibility.</p>
       <div class="scrollx"><table class="dt heat"><thead><tr><th>driver</th>%s</tr></thead><tbody>%s</tbody><tfoot><tr><th>driver</th>%s</tr></tfoot></table></div>',
      ghdr, grows, ghdr)
  }

  sprintf('
  <section class="branch" id="br-%s"%s>
    <div class="tiles">%s</div>
    %s
    %s
    %s
    %s
    %s
    %s
    %s
    %s
    %s
    %s
  </section>',
    sp, if (is_first) "" else " hidden",
    tiles, perf_html, ptable, re_heat, map_html, count_settings_html, comp_section, org_section, mnl_settings_html, gamma_section, trace_imgs)
}

sections <- paste(sapply(seq_along(species), function(i)
  build_species_section(species[i], i == 1)), collapse = "\n")

tabs <- paste(sprintf('<button class="tab%s" data-b="%s">%s</button>',
  ifelse(species == species[1], " active", ""), species, species), collapse = "")

# ---------------------------------------------------------------------------
# 9.  Model description + pipeline overview
# ---------------------------------------------------------------------------
model_desc <- '
<section class="model-desc">
<h2>Model architecture &amp; pipeline</h2>
<p class="note" style="border-left-color:var(--acc)">This report summarises a three-stage Bayesian livestock allocation model that converts EU-wide census data into spatially-explicit livestock densities at 10&thinsp;km resolution, disaggregated by species (BOV/SGT), production type (Dairy/Meat/Follower), and farming system (conventional/organic).</p>

<div class="pipeline">
<h3>Stage 1 &mdash; Totals count model <span class="sub">run_prior_module_count_model.R &rarr; run/count.R</span></h3>
<div class="formula-box">
<code>Y<sub>it</sub> ~ NegBin(&mu;<sub>it</sub>, r)</code><br>
<code>log(&mu;<sub>it</sub>) = <strong>X</strong><sub>it</sub> &middot; <strong>&beta;</strong><sub>g[i]</sub> + log(PastureArea<sub>it</sub>)</code><br>
<code>&beta;<sub>g</sub> ~ N(&mu;<sub>&beta;</sub>, &Sigma;<sub>RE</sub>)</code> &emsp; <em>(country random effects on intercept + key drivers)</em><br>
<code>&mu;<sub>&beta;</sub> ~ HS(&tau;<sub>0</sub>)</code> &emsp; <em>(horseshoe prior for automatic covariate selection)</em>
</div>
<table class="dt desc-table">
<tr><td class="cl">Outcome Y</td><td>LSU-equivalent livestock counts per NUTS3 region &times; year</td></tr>
<tr><td class="cl">Offset</td><td>log(total grazing area) &mdash; grazing-only exposure, summed over the Pasture/Grassland classes of whichever classification was fitted (the run log names them)</td></tr>
<tr><td class="cl">Grazing covariates</td><td><strong>Proportion method</strong>: the grazing-class shares, with the largest grazing class by area dropped as the reference</td></tr>
<tr><td class="cl">Non-pasture LU</td><td>log1p(area) for Cropland (6 classes), Forests (6), Natural_unmanaged, Urban</td></tr>
<tr><td class="cl">Biophysical</td><td>GDD5, Precipitation, Slope/Elevation SD, terrain shares (flat/steep/lowland/upland)</td></tr>
<tr><td class="cl">Socioeconomic</td><td>GHM_HI, GHM_TI, log1p(GDP), log1p(Pop), log1p(allPA_area)</td></tr>
<tr><td class="cl">Panel structure</td><td>3 years (2000, 2010, 2020) with year fixed effects (ref = 2000)</td></tr>
<tr><td class="cl">Random effects</td><td>Country-level on intercept + key socioeconomic &amp; terrain drivers &mdash; the grouping key is the column heading the heatmaps below</td></tr>
<tr><td class="cl">Input data</td><td><code>prior_model_1km_master_inputs.parquet</code> &rarr; aggregated to NUTS3</td></tr>
<tr><td class="cl">Script</td><td><code>run/count.R</code> &rarr; <code>run_prior_module_count_model.R</code></td></tr>
</table>

<h3>Stage 2 &mdash; D/O/F composition &delta; <span class="sub">composition/fit_composition.R</span></h3>
<div class="formula-box">
<code>(n<sub>D</sub>, n<sub>O</sub>, n<sub>F</sub>)<sub>j</sub> ~ MNL(softmax(<strong>X</strong><sub>j</sub> &middot; <strong>&delta;</strong><sub>k</sub>))</code> &emsp; <em>sum-to-zero symmetric MNL</em><br>
<code>&delta;<sub>k</sub> ~ HS</code> &emsp; <em>horseshoe for sparse, visible allocation drivers</em>
</div>
<table class="dt desc-table">
<tr><td class="cl">Outcome</td><td>Eurostat NUTS2 head counts: Dairy cows (D), Suckler/Meat (O), Followers (F)</td></tr>
<tr><td class="cl">Drivers</td><td>Same LU areas + biophysical (decorrelated, |cor| &lt; 0.70); LU classes protected from pruning</td></tr>
<tr><td class="cl">RE</td><td>Country-level on intercept</td></tr>
<tr><td class="cl">Data source</td><td><code>agr_r_animal</code> (Eurostat) + count-model dat_admin</td></tr>
<tr><td class="cl">Script</td><td><code>prep/prepare_composition_training.R</code> &rarr; <code>composition/fit_composition.R</code></td></tr>
</table>

<h3>Stage 3 &mdash; Organic/conventional split &delta;<sub>org</sub> <span class="sub">composition/fit_organic.R</span></h3>
<div class="formula-box">
<code>(n<sub>org</sub>, n<sub>conv</sub>)<sub>j</sub> ~ MNL(softmax(<strong>X</strong><sub>j</sub> &middot; &delta;<sub>org</sub>))</code><br>
<code>&delta;<sub>org</sub> ~ HS</code> &emsp; <em>spatial pattern from the organic-grassland proxy map</em><br>
<code>National level anchored to Eurostat BOV_org / SGT_org shares via per-country intercept c</code>
</div>
<table class="dt desc-table">
<tr><td class="cl">Outcome</td><td>Organic grassland proxy (Pasture_HIO + LIO vs HI + LI) at NUTS3</td></tr>
<tr><td class="cl">Anchor</td><td>Eurostat national organic livestock shares &rarr; per-country intercept offset c</td></tr>
<tr><td class="cl">Script</td><td><code>composition/fit_organic.R</code> &rarr; <code>composition/build_subclass_parameters.R</code></td></tr>
</table>

<h3>Final assembly &mdash; 12-cell &gamma; <span class="sub">composition/build_subclass_parameters.R</span></h3>
<div class="formula-box">
<code>&gamma;<sub>species,DOF,system,driver</sub> = &beta;<sub>total,country</sub> + &delta;<sub>DOF</sub> + &delta;<sub>org,signed</sub></code><br>
<code>cell weight<sub>i</sub> &prop; exp(&sum;<sub>d</sub> &gamma;<sub>d</sub> &middot; X<sub>id</sub> + offset<sub>i</sub>)</code> &emsp; <em>renormalized within NUTS3</em>
</div>
<table class="dt desc-table">
<tr><td class="cl">Cells</td><td>12 per species: {Dairy, Meat, Follower} &times; {conventional, organic}</td></tr>
<tr><td class="cl">Output</td><td><code>subclass_country_parameters.csv</code> &mdash; per (country, subclass, system, driver)</td></tr>
<tr><td class="cl">Downscale</td><td>Predict on 10 km &times; NUTS3 grid &rarr; reconcile to NUTS3 admin totals</td></tr>
</table>

<h3>Shrinkage &amp; Regularization Architecture</h3>
<p class="note" style="border-left-color:#d97706;background:#fffbeb">To prevent overfitting across dozens of spatial drivers while maintaining model expressiveness, Bayesian shrinkage is applied systematically across all tiers of the model hierarchy.</p>
<table class="dt desc-table">
<thead>
<tr>
  <th>Model Tier</th>
  <th>Shrinkage Mechanism</th>
  <th>Target Parameters (Where Applied)</th>
  <th>Unshrunk Parameters</th>
  <th>Hyperparameters &amp; Behavior</th>
</tr>
</thead>
<tbody>
<tr>
  <td class="cl"><strong>Count Model Slopes</strong></td>
  <td><span class="badge" style="background:#e0f2fe;color:#0369a1">Regularized Horseshoe</span><br>(Piironen &amp; Vehtari 2017)</td>
  <td>All 36 continuous drivers: climate, terrain, pasture proportions, non-pasture LU areas, year FEs</td>
  <td><span class="badge" style="background:#f3f4f6;color:#374151">Global Intercept</span><br>(flat diffuse prior)</td>
  <td>&tau;<sub>0</sub> = 0.50, slab c<sup>2</sup> = 100, local &lambda;<sub>j</sub> ~ C<sup>+</sup>(0,1). Shrinkage factor &kappa;<sub>j</sub> &in; [0,1] pulls noise drivers strictly to 0 while allowing genuine signals past the Finnish slab cap.</td>
</tr>
<tr>
  <td class="cl"><strong>Country Gatekeeper &amp; RE Variance</strong></td>
  <td><span class="badge" style="background:#fef3c7;color:#92400e">Hierarchical Gatekeeper &tau;<sub>m</sub> + Regularized Slab</span></td>
  <td>Country-level random slope deviations <strong>u</strong><sub>m</sub> across 35 drivers</td>
  <td><span class="badge" style="background:#f3f4f6;color:#374151">Country Intercepts &mu;<sub>0</sub> + u<sub>m,0</sub></span><br>(independent baseline variation)</td>
  <td>&tau;<sub>m</sub> ~ C<sup>+</sup>(0, 1.0) gates slope deviations per country. Low-observation countries (LU, MT, CY) automatically shrink &tau;<sub>m</sub> &rarr; 0 (forcing full pooling to continental mean &mu;<sub>&beta;</sub>), while high-observation countries (DE, FR, ES) retain flexibility &tau;<sub>m</sub> &approx; 1. Base variance capped at &sigma;<sub>RE</sub> &le; 10 (c<sup>2</sup>=100).</td>
</tr>
<tr>
  <td class="cl"><strong>Subtype Composition &delta;<sub>DOF</sub></strong></td>
  <td><span class="badge" style="background:#f3e8ff;color:#7e22ce">Symmetric Horseshoe</span></td>
  <td>27 allocation drivers across Dairy (D), Meat (O), and Follower (F) categories</td>
  <td>Subclass country intercepts</td>
  <td>Coupled cross-equation horseshoe prior (M<sub>sym</sub>) with sum-to-zero constraint &sum; &delta;<sub>k</sub> = 0 and estimated slab c<sup>2</sup> ~ InvGamma(20,4).</td>
</tr>
<tr>
  <td class="cl"><strong>Organic Contrast &delta;<sub>org</sub></strong></td>
  <td><span class="badge" style="background:#dcfce7;color:#15803d">Symmetric Horseshoe + National Anchor</span></td>
  <td>23 spatial drivers of organic vs conventional farming systems</td>
  <td>National level anchor intercepts</td>
  <td>Regularized spatial gradient from organic grassland proxies, with country-level intercepts <em>c</em> exactly anchored to Eurostat national statistics.</td>
</tr>
</tbody>
</table>
</div>
</section>'

# ---------------------------------------------------------------------------
# 10. CSS
# ---------------------------------------------------------------------------
css <- '
:root{--bg:#FAFBFA;--surf:#fff;--ink:#14211E;--mut:#5E6E69;--line:#E5EAE7;--acc:#0F7A67;--acc2:#075447;--warn:#B26B00;--warnbg:#FBF3E4;--bad:#B2182B;--good:#0F7A67}
@media(prefers-color-scheme:dark){:root{--bg:#0E1613;--surf:#14201C;--ink:#E8F0EC;--mut:#93A39D;--line:#25332E;--acc:#54BCA6;--acc2:#8AD3C3;--warn:#E0A050;--warnbg:#241d10;--bad:#E27D8C}}
:root[data-theme=dark]{--bg:#0E1613;--surf:#14201C;--ink:#E8F0EC;--mut:#93A39D;--line:#25332E;--acc:#54BCA6;--acc2:#8AD3C3;--warn:#E0A050;--warnbg:#241d10;--bad:#E27D8C}
:root[data-theme=light]{--bg:#FAFBFA;--surf:#fff;--ink:#14211E;--mut:#5E6E69;--line:#E5EAE7;--acc:#0F7A67;--acc2:#075447;--warn:#B26B00;--warnbg:#FBF3E4;--bad:#B2182B;--good:#0F7A67}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font-family:ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;line-height:1.5;font-size:15px}
.wrap{max-width:1160px;margin:0 auto;padding:0 22px 80px}
header.top{position:sticky;top:0;z-index:20;background:color-mix(in srgb,var(--bg) 88%,transparent);backdrop-filter:blur(8px);border-bottom:1px solid var(--line)}
.hd{max-width:1160px;margin:0 auto;padding:14px 22px 0}
h1{font-size:22px;font-weight:640;letter-spacing:-.01em;margin:0 0 2px;text-wrap:balance}
h2{font-size:18px;font-weight:640;margin:32px 0 12px;letter-spacing:-.01em;border-bottom:2px solid var(--line);padding-bottom:6px}
.subttl{color:var(--mut);font-size:13.5px;margin:0 0 12px}
.tabs{display:flex;gap:4px;flex-wrap:wrap}
.tab{appearance:none;border:1px solid var(--line);border-bottom:none;background:transparent;color:var(--mut);padding:9px 15px 10px;border-radius:9px 9px 0 0;cursor:pointer;font:inherit;font-weight:600}
.tab.active{color:var(--ink);background:var(--surf);border-color:var(--line);box-shadow:0 -2px 0 var(--acc) inset}
.branch{padding-top:20px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:6px 0 20px}
.tile{background:var(--surf);border:1px solid var(--line);border-radius:12px;padding:13px 14px}
.tval{font-size:23px;font-weight:650;letter-spacing:-.02em;font-variant-numeric:tabular-nums}
.tlabel{color:var(--mut);font-size:12px;margin-top:2px}
.tsub{color:var(--mut);font-size:11px;opacity:.8;margin-top:3px}
h3{font-size:15.5px;font-weight:620;margin:30px 0 10px;letter-spacing:-.005em}
h3 .sub{font-weight:400;font-size:12px;color:var(--mut);margin-left:8px}
.note{color:var(--mut);font-size:12.7px;background:var(--surf);border:1px solid var(--line);border-left:3px solid var(--acc);border-radius:8px;padding:10px 13px;margin:0 0 12px}
table.dt{width:100%;border-collapse:collapse;font-size:12.8px}
table.dt th,table.dt td{padding:5px 9px;border-bottom:1px solid var(--line);text-align:left}
table.dt th{color:var(--mut);font-weight:600;font-size:11px;text-transform:uppercase;letter-spacing:.03em}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
td.cl{font-weight:500}
td.ci{color:var(--mut);font-size:11.5px}
td.sig{text-align:center;font-size:14px}
tr.pos td{background:color-mix(in srgb,var(--good) 8%,transparent)} tr.pos td.sig{color:var(--good)}
tr.neg td{background:color-mix(in srgb,var(--bad) 8%,transparent)}  tr.neg td.sig{color:var(--bad)}
tr.ns td.sig{color:var(--mut);opacity:.3}
.badge{display:inline-block;font-size:9px;font-weight:700;text-transform:uppercase;letter-spacing:.04em;padding:1px 5px;border-radius:4px;margin-left:6px;vertical-align:middle}
.badge.prop{background:#E0F5F0;color:#0F7A67} .badge.lu{background:#F0EBF8;color:#6B4EA0} .badge.yr{background:#FBF3E4;color:#B26B00}
@media(prefers-color-scheme:dark){.badge.prop{background:#14302A;color:#54BCA6}.badge.lu{background:#24193A;color:#B29FD6}.badge.yr{background:#241d10;color:#E0A050}}
:root[data-theme=dark] .badge.prop{background:#14302A;color:#54BCA6}:root[data-theme=dark] .badge.lu{background:#24193A;color:#B29FD6}:root[data-theme=dark] .badge.yr{background:#241d10;color:#E0A050}
.scrollx{overflow-x:auto;border:1px solid var(--line);border-radius:10px;background:#fff}
table.heat{font-size:11px;border-collapse:separate;border-spacing:0;background:#fff;color:#14211E}
table.heat th.rot{height:96px;white-space:nowrap;vertical-align:bottom;padding:0}
table.heat th.rot span{display:inline-block;transform:rotate(-90deg);transform-origin:left;translate:12px -6px;font-size:10px;color:#5E6E69}
table.heat td.hm{position:relative;padding:0;height:24px;vertical-align:middle;text-align:center;font-variant-numeric:tabular-nums;color:#0a1512;border-bottom:1px solid rgba(255,255,255,.5)}
table.heat .tbg{position:absolute;top:0;left:0;right:0;bottom:0;background:#f9f9f9;z-index:0}
table.heat .tile{position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);z-index:1;border-radius:1px}
table.heat .hval{position:relative;z-index:2;display:block;padding:0 4px;font-size:10.5px}
table.heat td.sticky{position:sticky;left:0;background:#fff;color:#14211E;font-weight:600}
figure{margin:0}
figure.mapc img{width:100%;height:auto;display:block;border:1px solid var(--line);border-radius:10px;background:var(--surf)}
.figure{margin:16px 0 24px}
.figure h4{margin:0 0 6px;font-size:14px;font-weight:600}
.mapgrid{display:grid;grid-template-columns:repeat(auto-fill,minmax(330px,1fr));gap:12px}
footer{color:var(--mut);font-size:12px;margin-top:40px;border-top:1px solid var(--line);padding-top:16px}
.themebtn{position:fixed;right:14px;bottom:14px;z-index:30;border:1px solid var(--line);background:var(--surf);color:var(--ink);border-radius:999px;padding:8px 12px;cursor:pointer;font:inherit;font-size:12px}
.formula-box{background:var(--surf);border:1px solid var(--line);border-radius:10px;padding:12px 16px;margin:8px 0 14px;font-family:"SF Mono",Menlo,Consolas,monospace;font-size:13px;line-height:1.8;overflow-x:auto}
.formula-box code{font-size:13px}
table.desc-table td.cl{width:170px;color:var(--mut);font-size:12px;text-transform:uppercase;letter-spacing:.02em}
table.desc-table td{font-size:13px}
table.desc-table code{font-size:12px;background:color-mix(in srgb,var(--acc) 10%,transparent);padding:1px 4px;border-radius:3px}
.pipeline h3{margin-top:22px}
.model-desc{margin-bottom:10px}
'

js <- "
document.querySelectorAll('.tab').forEach(t=>t.addEventListener('click',()=>{
 document.querySelectorAll('.tab').forEach(x=>x.classList.remove('active'));t.classList.add('active');
 document.querySelectorAll('.branch').forEach(s=>s.hidden=true);
 document.getElementById('br-'+t.dataset.b).hidden=false;window.scrollTo({top:0,behavior:'smooth'});}));
const tb=document.getElementById('themebtn');function cur(){return document.documentElement.getAttribute('data-theme')|| (matchMedia('(prefers-color-scheme:dark)').matches?'dark':'light')}
tb.addEventListener('click',()=>{const n=cur()==='dark'?'light':'dark';document.documentElement.setAttribute('data-theme',n);tb.textContent=n==='dark'?'\\u2600 light':'\\u263E dark'});
"

html <- sprintf('<!doctype html><html><head><meta charset="utf-8">
<title>Livestock count model &mdash; fit report (%s)</title>
<style>%s</style></head><body>
<header class="top"><div class="hd"><h1>Livestock count model &mdash; parameter &amp; fit report (%s)</h1>
<p class="subttl">NUTS3 admin level &middot; Negative-binomial with country RE &middot; GLOBIOM classification &middot; Pasture proportion method (ref: Pasture_HI)</p>
<div class="tabs">%s</div></div></header>
<div class="wrap">
%s
%s
<footer>Generated %s &middot; live production posterior means &middot; pasture-only offset (log total pasture area) &middot; composition &#948; from NUTS2 Eurostat head counts</footer>
</div>
<button class="themebtn" id="themebtn">&#9790; dark</button>
<script>%s</script>
</body></html>', format(Sys.Date(), "%Y-%m-%d"), css, format(Sys.Date(), "%Y-%m-%d"), tabs, model_desc, sections, format(Sys.Date()), js)

writeLines(html, OUTHTML)
cat("wrote", OUTHTML, "\n")
cat("size:", round(file.info(OUTHTML)$size / 1e6, 1), "MB\n")

# Also copy to count_model_fit_report.html
file.copy(OUTHTML, "output/report/count_model_fit_report.html", overwrite = TRUE)
cat("copied to output/report/count_model_fit_report.html\n")
