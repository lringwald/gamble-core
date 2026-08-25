#!/usr/bin/env Rscript
# =============================================================================
# nested_report.R — self-contained HTML overview of a NESTED fit, ONE FILE PER FIT
# =============================================================================
# Usage:  Rscript postprocess/nested_report.R [fit1.rds fit2.rds ...]
#         (no arguments -> the newest output/nested_cut_*.rds)
# Writes: output/report/nested_fit_<fit name>.html   (one per fit, images inlined)
#
# Deliberately NO cross-fit comparison: each report describes a single fit on its own terms.
#
# WHY A SEPARATE ENGINE. postprocess/engine.R reconstructs the FLAT model from on-disk posterior
# batches + model_metadata.qs; a nested fit has neither -- it is a tree of per-node draws. This
# reads the saved fit and predicts through the tree (predict_nested_cut).
#
# THE SCORES ARE IN-SAMPLE: reproduction, not generalisation. Held-out comparison lives in
# experiments/mixing/re_idx_tradeoff.R.
# =============================================================================
suppressMessages({library(qs2); library(data.table); library(ggplot2); library(scales); library(base64enc)})
source("codes/mnlogit_rcpp_sym.R"); source("codes/nested_cut.R"); source("codes/nest_trees.R")
`%||%` <- function(a, b) if (is.null(a)) b else a

FITS <- commandArgs(trailingOnly = TRUE)
if (!length(FITS)) FITS <- head(sort(Sys.glob("output/nested_cut_*.rds"), decreasing = TRUE), 1)
FITS <- FITS[file.exists(FITS)]
if (!length(FITS)) stop("no nested fits found")
N_MAP_CLASSES <- as.integer(Sys.getenv("REPORT_MAP_CLASSES", "12"))
dir.create("output/report", showWarnings = FALSE, recursive = TRUE)
MAPDIR <- "output/report_work/maps_nested"; dir.create(MAPDIR, showWarnings = FALSE, recursive = TRUE)

# reuse the flat report's stylesheet so both reports look like one product
CSS <- local({ e <- new.env()
  for (ex in parse("postprocess/build_html.R"))
    if (is.call(ex) && identical(as.character(ex[[1]]), "<-") &&
        is.name(ex[[2]]) && identical(as.character(ex[[2]]), "css")) eval(ex, e)
  e$css %||% "body{font-family:system-ui;margin:2rem}" })
EXTRA_CSS <- paste0("\n/* supplement: rules the shared stylesheet does not define */\n",
  "figure{margin:14px 0}figure img{width:100%;border:1px solid var(--line);border-radius:6px}\n",
  "section h3{margin:18px 0 6px;font-size:14px;color:var(--acc2)}\n",
  "div.bar{background:var(--line);border-radius:3px;height:8px;width:120px}\n",
  "div.bar i{display:block;height:8px;border-radius:3px}\n")

esc  <- function(s) { s <- gsub("&", "&amp;", s); s <- gsub("<", "&lt;", s); gsub(">", "&gt;", s) }
fnum <- function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = ",")
tile <- function(label, val, sub = "") sprintf(
  '<div class="tile"><div class="tval">%s</div><div class="tlabel">%s</div>%s</div>',
  val, esc(label), if (nzchar(sub)) sprintf('<div class="tsub">%s</div>', esc(sub)) else "")
b64 <- function(p) if (!file.exists(p)) "" else
  paste0("data:image/png;base64,", base64enc::base64encode(readBin(p, "raw", file.info(p)$size)))

seqcols <- c("#F2FbF8","#CDEBE3","#98D6C8","#5FBBA8","#2E9A86","#0F7A67","#075447")
theme_map <- theme_void(base_size = 11) + theme(
  plot.title = element_text(size = 11, hjust = .5, margin = margin(b = 3)),
  legend.position = "right", legend.key.height = unit(14, "pt"), legend.key.width = unit(8, "pt"),
  plot.background = element_rect(fill = "white", colour = NA),
  strip.text = element_text(size = 10, face = "bold"))

score_fit <- function(f) {
  z <- readRDS(f); inp <- readRDS(z$inputs)
  X <- as.matrix(inp$X_mat); if (is.null(colnames(X))) colnames(X) <- inp$col_names
  X[!is.finite(X)] <- 0
  Y <- as.matrix(inp$Y_pixel); Y <- Y / rowSums(Y); Y <- Y[, z$cats, drop = FALSE]
  # key on the INTEGER group_idx: the fit's group_levels are unique(group_idx) in APPEARANCE ORDER,
  # while z$group holds country NAMES -- passing names silently degrades to pooled prediction.
  P <- predict_nested_cut(z$fit, X, D = z$fit$M, summary = "mean", group_idx = inp$group_idx_vec)
  ll <- sum(Y * log(pmax(P, 1e-12)))
  s0 <- colSums(Y) / sum(Y)
  ll0 <- sum(Y * log(matrix(s0, nrow(Y), ncol(Y), byrow = TRUE)))
  cn <- inp$class_nest
  nst <- if (is.null(cn)) setNames(rep("(flat)", length(z$cats)), z$cats)
         else setNames(as.character(cn[[2]])[match(z$cats, cn[[1]])], z$cats)
  message(sprintf("    %-52s McF %.4f", basename(f), 1 - ll / ll0))
  list(file = f, z = z, Y = Y, P = P, ll = ll, ll0 = ll0, mcf = 1 - ll / ll0,
       mad = mean(abs(P - Y)), coordX = z$coordX, coordY = z$coordY, nest = nst,
       comp = data.table(class = z$cats, obs = 100 * s0, prd = 100 * colSums(P) / sum(P),
                         ratio = (colSums(P) / sum(P)) / s0))
}

render_maps <- function(r) {
  tops <- r$comp[order(-obs)][seq_len(min(N_MAP_CLASSES, .N)), class]
  tag <- tools::file_path_sans_ext(basename(r$file)); out <- setNames(character(0), character(0))
  for (cl in tops) {
    d <- rbind(data.table(x = r$coordX, y = r$coordY, v = r$Y[, cl], panel = "observed"),
               data.table(x = r$coordX, y = r$coordY, v = r$P[, cl], panel = "predicted"))
    pl <- ggplot(d, aes(x, y, fill = v)) + geom_raster() + facet_wrap(~panel, nrow = 1) +
      scale_fill_gradientn(colours = seqcols, labels = percent_format(accuracy = 1), name = NULL,
                           limits = c(0, max(d$v, na.rm = TRUE))) +
      coord_equal() + labs(title = cl) + theme_map
    f <- file.path(MAPDIR, sprintf("%s__%s.png", tag, gsub("[^A-Za-z0-9_]", "_", cl)))
    suppressWarnings(ggsave(f, pl, width = 7.2, height = 3.4, dpi = 96))
    out[cl] <- f
  }
  out
}

body_html <- function(r) {
  z <- r$z
  cp <- copy(r$comp)[order(-obs)][, nest := r$nest[class]]
  bar <- function(x) { w <- max(0, min(100, 50 * x))
    col <- if (abs(log(x)) < .15) "var(--acc)" else if (abs(log(x)) < .4) "var(--warn)" else "var(--bad)"
    sprintf('<div class="bar"><i style="width:%.0f%%;background:%s"></i></div>', w, col) }
  rows <- paste0(vapply(seq_len(nrow(cp)), function(i) sprintf(
    '<tr><td class="cl">%s</td><td class="ci">%s</td><td class="n">%.3f</td><td class="n">%.3f</td><td class="n">%.2f</td><td>%s</td></tr>',
    esc(cp$class[i]), esc(cp$nest[i] %||% ""), cp$obs[i], cp$prd[i], cp$ratio[i], bar(cp$ratio[i])), character(1)), collapse = "")
  nt <- cp[, .(obs = sum(obs), prd = sum(prd)), by = nest][order(-obs)][, ratio := prd / obs]
  nrows <- paste0(vapply(seq_len(nrow(nt)), function(i) sprintf(
    '<tr><td class="cl">%s</td><td class="n">%.2f</td><td class="n">%.2f</td><td class="n">%.3f</td></tr>',
    esc(nt$nest[i]), nt$obs[i], nt$prd[i], nt$ratio[i]), character(1)), collapse = "")
  prov <- sprintf('<table class="dt"><tbody>%s</tbody></table>', paste0(sprintf(
    '<tr><th>%s</th><td class="cl">%s</td></tr>',
    c("design", "design written", "tree", "variant", "RE block", "shrinkage", "classes", "chains"),
    c(esc(basename(z$inputs %||% "?")), esc(as.character(z$inputs_mtime %||% "?")), esc(z$tree_source %||% "?"),
      if (isTRUE(z$use_iv)) "iv-nested (lambda estimated)" else "factorized (levels independent)",
      esc(z$re_cols %||% "?"),
      if (isTRUE(z$symmetric_hs)) "symmetric (zero-sum/CLR) horseshoe" else "diagonal horseshoe",
      length(z$cats), z$n_chains %||% 1)), collapse = ""))
  lam <- ""
  if (isTRUE(z$use_iv)) { idt <- tryCatch(nested_cut_identification(z$fit), error = function(e) NULL)
    if (!is.null(idt)) lam <- paste0('<h3>Inclusive-value &lambda;</h3><table class="dt"><thead><tr>',
      '<th>node</th><th>child</th><th class="n">IV~design R2</th><th class="n">&lambda;</th><th class="n">95% CI</th></tr></thead><tbody>',
      paste0(sprintf('<tr><td class="cl">%s</td><td>%s</td><td class="n">%.2f</td><td class="n">%.3f</td><td class="n">[%.3f, %.3f]</td></tr>',
        idt$node, idt$iv_child, idt$iv_r2, idt$lambda, idt$lambda_q025, idt$lambda_q975), collapse = ""),
      '</tbody></table><p class="note">&lambda; outside (0,1] has no random-utility interpretation; R2 near 1 means &lambda; is not separately identified from the direct effects.</p>') }
  sprintf('<section><h2>Fit overview</h2><div class="tiles">%s%s%s%s</div>%s%s
<h3>Composition by nest</h3><table class="dt"><thead><tr><th>nest</th><th class="n">observed %%</th><th class="n">predicted %%</th><th class="n">ratio</th></tr></thead><tbody>%s</tbody></table>
<h3>Composition by class</h3><table class="dt"><thead><tr><th>class</th><th>nest</th><th class="n">obs %%</th><th class="n">pred %%</th><th class="n">ratio</th><th></th></tr></thead><tbody>%s</tbody></table></section>',
    tile("McFadden", sprintf("%.4f", r$mcf), "in-sample"),
    tile("log-likelihood", fnum(r$ll, 0), sprintf("null %s", fnum(r$ll0, 0))),
    tile("mean abs. error", sprintf("%.5f", r$mad), "share units"),
    tile("classes", length(z$cats), if (isTRUE(z$symmetric_hs)) "symmetric HS" else "diagonal HS"),
    prov, lam, nrows, rows)
}

write_report <- function(r) {
  mp <- render_maps(r)
  maps <- paste0('<section><h2>Observed vs predicted</h2><p class="note">Top ', length(mp),
    ' classes by area. Left: observed share, right: posterior-mean prediction; shared colour scale within each class.</p>',
    paste0(vapply(names(mp), function(cl) sprintf('<figure><img src="%s" alt="%s"></figure>',
      b64(mp[cl]), esc(cl)), character(1)), collapse = ""), '</section>')
  caveat <- paste0('<div class="caveat"><b>These scores are in-sample.</b> This fit trained on all rows, ',
    'so they measure how well it reproduces its own data, not how it generalises. For comparison between ',
    'configurations use <code>experiments/mixing/re_idx_tradeoff.R</code>, which holds out a country-stratified 25%.</div>')
  html <- sprintf('<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Nested fit - %s</title>
<style>%s</style></head><body>
<header class="top"><div class="hd"><h1>Nested land-use model - fit report</h1>
<p class="subttl">%s</p></div></header>
<div class="wrap">%s%s%s
<footer>Predicted through the tree with <code>predict_nested_cut</code> (per-group random effects, %d draws), generated %s.</footer>
</div></body></html>',
    esc(basename(r$file)), paste0(CSS, EXTRA_CSS), esc(basename(r$file)),
    caveat, body_html(r), maps, r$z$fit$M, format(Sys.time(), "%Y-%m-%d %H:%M"))
  f <- file.path("output/report", paste0("nested_fit_", tools::file_path_sans_ext(basename(r$file)), ".html"))
  writeLines(html, f)
  cat(sprintf("wrote %s (%.1f MB)\n", f, file.info(f)$size / 1e6))
}

message(">>> ", length(FITS), " fit(s)")
for (f in FITS) write_report(score_fit(f))
