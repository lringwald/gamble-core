# =============================================================================
# HTML validation report for the livestock count model.
# Reads output/count_validation.rds (produced by DRIVER_VALIDATE=TRUE) + the totals convergence checks
# + the 12-cell subclass parameters, and writes a self-contained HTML.
#   Rscript postprocess/count_validation_report.R
# =============================================================================
suppressMessages(library(data.table))
VF <- "output/count_validation.rds"
if (!file.exists(VF)) stop("No ", VF, " — run the driver with DRIVER_VALIDATE=TRUE first.")
v <- readRDS(VF)
OUT <- "output/report/count_validation_report.html"; dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

fmt <- function(x, d = 3) formatC(x, digits = d, format = "f")
rows_metrics <- ""; rows_year <- ""; rows_ds <- ""
for (sp in names(v)) {
  m <- v[[sp]]$metrics
  rows_metrics <- paste0(rows_metrics, sprintf(
    "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td class='oos'>%s</td><td class='oos'>%s</td></tr>",
    sp, fmt(m["in_sample","cor_log"]), fmt(m["in_sample","R2_log"]), fmt(m["in_sample","cor"]),
    fmt(m["oos","cor_log"]), fmt(m["oos","R2_log"])))
  # by-year in/out log-cor
  Y <- v[[sp]]$Y; yr <- v[[sp]]$year
  for (yy in sort(unique(yr))) { r <- yr == yy & is.finite(Y)
    icl <- suppressWarnings(cor(log1p(pmax(v[[sp]]$mu_in[r],0)), log1p(Y[r])))
    ocl <- suppressWarnings(cor(log1p(pmax(v[[sp]]$mu_oos[r],0)), log1p(Y[r])))
    rows_year <- paste0(rows_year, sprintf("<tr><td>%s</td><td>%d</td><td>%s</td><td class='oos'>%s</td><td>%d</td></tr>",
      sp, yy, fmt(icl), fmt(ocl), sum(r))) }
  if (!is.null(v[[sp]]$downscale)) { d <- v[[sp]]$downscale
    agg <- d[, .(grid_total = sum(mu, na.rm = TRUE)), by = year]
    for (i in seq_len(nrow(agg))) rows_ds <- paste0(rows_ds, sprintf(
      "<tr><td>%s</td><td>%d</td><td>%s units</td><td>%.0f</td></tr>", sp, agg$year[i], formatC(nrow(d[year==agg$year[i]]), big.mark=","), agg$grid_total[i]))
  }
}
# totals convergence (max Rhat) if the check files are present
conv <- ""
for (f in Sys.glob("output/count_rcpp_admin_convergence_check_*.rds")) {
  cv <- tryCatch(as.data.table(readRDS(f)), error = function(e) NULL); if (is.null(cv)) next
  rc <- grep("rhat", names(cv), ignore.case = TRUE, value = TRUE)[1]; if (is.na(rc)) next
  rh <- cv[[rc]]; rh <- rh[is.finite(rh)]; sp <- sub(".*production_([A-Z]+)_.*", "\\1", basename(f))
  conv <- paste0(conv, sprintf("<tr><td>%s</td><td>%s</td><td>%.0f%%</td><td>%.0f%%</td></tr>",
    sp, fmt(max(rh)), 100*mean(rh<1.01), 100*mean(rh<1.1))) }

css <- "body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;max-width:1000px;margin:2rem auto;padding:0 1rem;color:#1a2230;line-height:1.5}
h1{font-size:1.6rem}h2{font-size:1.15rem;margin-top:2rem;border-bottom:2px solid #e3e8f0;padding-bottom:.3rem}
table{border-collapse:collapse;width:100%;margin:.6rem 0;font-variant-numeric:tabular-nums}
th,td{border:1px solid #dde3ec;padding:.4rem .6rem;text-align:right}th:first-child,td:first-child{text-align:left}
th{background:#f2f5fa}.oos{color:#8a5a00;background:#fff8ec}.note{color:#5b6472;font-size:.9rem}"
html <- paste0("<!doctype html><html><head><meta charset='utf-8'><title>Count model validation</title><style>",css,
"</style></head><body><h1>Livestock count model — validation</h1>",
"<p class='note'>In-sample uses the country random effect; the <span class='oos'>OOS proxy</span> is the pooled (no-RE) prediction — an unseen-country approximation. Metrics on a log scale (counts are heavy-tailed). Offset = log(area) exposure.</p>",
"<h2>Totals convergence</h2><table><tr><th>species</th><th>max Rhat</th><th>%&lt;1.01</th><th>%&lt;1.1</th></tr>", conv, "</table>",
"<h2>Fit: in-sample (RE) vs OOS proxy (pooled)</h2><table><tr><th>species</th><th>cor(log)</th><th>R²(log)</th><th>cor</th><th class='oos'>OOS cor(log)</th><th class='oos'>OOS R²(log)</th></tr>", rows_metrics, "</table>",
"<h2>By year</h2><table><tr><th>species</th><th>year</th><th>in cor(log)</th><th class='oos'>OOS cor(log)</th><th>n</th></tr>", rows_year, "</table>",
"<h2>Downscale (10km × NUTS3 grid)</h2><table><tr><th>species</th><th>year</th><th>grid units</th><th>predicted total (×LSU unit)</th></tr>", rows_ds, "</table>",
"<p class='note'>Downscale predicts density per grid unit and aggregates to NUTS3; reconcile to the NUTS3 totals and split by the 12-cell subclass γ for the final allocation. A <b>pasture+cropland-area offset</b> (DRIVER_OFFSET=agri) allocates animals to grazing land more directly than total area.</p>",
"</body></html>")
writeLines(html, OUT)
cat(">>> wrote", OUT, "\n")
