#!/usr/bin/env Rscript
# =============================================================================
# run_all.R — one entry point for every test in the repo, and a report of it
# =============================================================================
# Why this exists: there was no single command that ran everything, so the pixel suite went
# unverified for a full day of sampler changes while targeted tests passed. Each file below is
# self-contained and exits non-zero on failure; this runs them, prints one table, and writes a
# self-contained HTML report with a figure.
#
#   Rscript tests/run_all.R              # everything + report
#   Rscript tests/run_all.R fast         # skip the two long suites
#   Rscript tests/run_all.R no-report    # table only, no PNG/HTML
#
# THREE STATUSES, not two:
#   PASS  assertions ran and all held
#   FAIL  a non-zero exit, a [FAIL] line, or NO ASSERTIONS AT ALL (see below)
#   SKIP  the file declared `SKIP:` on stdout because its inputs are absent
#
# A file that asserts NOTHING is a FAIL, not a pass. Exit code 0 with zero assertions means it did
# not run (unreadable file, early return, a stub) -- silent success is the failure mode this repo
# keeps hitting, and a green runner that hides it is worse than no runner. SKIP exists so a test
# whose fixtures are legitimately absent can say so OUT LOUD instead of being mistaken for either.
#
# THE RUN LIST IS ORDER ONLY. It used to be the whole membership list, hand-maintained, and three
# feature tests (re_mean_shift, joint_shrink, mundlak) sat outside it and never ran. Now every
# tests/test_*.R is discovered and RUN regardless; the list only fixes the order of the ones that
# should go first, and anything undeclared is reported so the omission is visible rather than silent.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
FAST      <- "fast" %in% args
NO_REPORT <- "no-report" %in% args
root <- if (basename(getwd()) == "tests") ".." else "."
setwd(root)

# Order-only preference; membership is by discovery (see header).
preferred <- c(
  "tests/test_nested_cut_units.R",    # fast helper units (M=1 IV fields) -- runs first, fails fast
  "tests/test_altspec_sampler.R",     # delta: recovery, multi-block, per_class, symmetric, scaling
  "tests/test_altspec_nested.R",      # blocks through nested_cut: wiring, IV path, interactions
  "tests/test_nested_cut.R",
  "tests/test_nested_cut_focal.R",
  "tests/test_nested_iv.R",
  "tests/test_re_mean_shift.R",       # sum-to-zero RE identification (mu == population average)
  "tests/test_joint_shrink.R",        # joint FE/RE gate kappa_v, incl. the seed-noise floor
  "tests/test_mundlak.R",             # group-mean prior mean: gamma recovery + beta debiasing
  "tests/test_predict_gamble.R",      # unified predictor: dispatch, delta, mundlak, disk round-trip
  "tests/test_target_class_split.R"   # project target cascade (SKIPs unless the BMLEH build exists)
)
long_suites <- c("tests/test_suite_lu_pixel.R", "tests/test_suite_ls_count.R")
# MEMBERSHIP IS BY DISCOVERY (see header). This line was lost in the 2026-09-17 restructure, leaving
# `discovered` assigned from itself: the file still PARSED, so the breakage surfaced only at runtime
# as "object 'discovered' not found" -- i.e. the whole runner was dead while looking healthy.
discovered <- sort(list.files("tests", pattern = "^test_.*\\.R$", full.names = TRUE))
discovered <- setdiff(discovered, long_suites)   # the long gates are appended separately, below
undeclared <- setdiff(discovered, preferred)
if (length(undeclared))
  cat(sprintf("  [note] running %d test file(s) not in the ordered list: %s\n",
              length(undeclared), paste(basename(undeclared), collapse = ", ")))
files <- c(preferred[preferred %in% discovered], undeclared)

# The two long gates live in tests/ and are self-contained (they locate SOURCE files, not fixtures),
# so they always run unless FAST is requested -- they never SKIP.
if (!FAST) files <- c(files, long_suites)
files <- files[file.exists(files)]

cat(sprintf("\n%s\n running %d test file(s)%s\n%s\n", strrep("=", 64), length(files),
            if (FAST) "  [fast: long suites skipped]" else "", strrep("=", 64)))
res <- data.frame(file = character(), status = character(), secs = numeric(),
                  npass = integer(), nfail = integer(), note = character(), stringsAsFactors = FALSE)
fail_lines <- list()
for (f in files) {
  t0 <- Sys.time()
  out <- suppressWarnings(system2("Rscript", f, stdout = TRUE, stderr = TRUE))
  rc <- attr(out, "status"); rc <- if (is.null(rc)) 0L else as.integer(rc)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # a suite reports "RESULT: a/b"; a targeted test reports [PASS]/[FAIL] lines
  hit    <- grep("RESULT:", out, value = TRUE)
  skipln <- grep("^\\s*SKIP:", out, value = TRUE)
  npass  <- sum(grepl("\\[PASS\\]", out)); nfail <- sum(grepl("\\[FAIL\\]", out))
  silent <- !length(hit) && npass == 0L && nfail == 0L && !length(skipln)
  status <- if (length(skipln) && rc == 0L && nfail == 0L) "SKIP"
            else if (rc != 0L || nfail > 0L || silent || (length(hit) && grepl("FAILED", hit[1]))) "FAIL"
            else "PASS"
  note <- if (status == "SKIP") trimws(sub(".*SKIP:", "", skipln[1]))
          else if (length(hit)) trimws(sub(".*RESULT:", "", hit[1]))
          else if (silent) "NO ASSERTIONS RUN -- did the file execute?"
          else sprintf("%d pass / %d fail", npass, nfail)
  res <- rbind(res, data.frame(file = basename(f), status = status, secs = el,
                               npass = npass, nfail = nfail, note = note, stringsAsFactors = FALSE))
  cat(sprintf("  [%s] %-34s %6.1fs  %s\n", status, basename(f), el, note))
  if (status == "FAIL") {
    fl <- tail(grep("\\[FAIL\\]|Error", out, value = TRUE), 3)
    fail_lines[[basename(f)]] <- fl
    if (length(fl)) cat(paste0("        ", fl, collapse = "\n"), "\n")
  }
}
nf <- sum(res$status == "FAIL"); ns <- sum(res$status == "SKIP")
cat(sprintf("\n%s\n %s  (%d file(s), %d assertions, %.1f min)%s\n%s\n", strrep("=", 64),
            if (nf == 0L) "ALL PASS" else sprintf("%d FILE(S) FAILED", nf),
            nrow(res), sum(res$npass) + sum(res$nfail), sum(res$secs) / 60,
            if (ns) sprintf("  [%d skipped]", ns) else "", strrep("=", 64)))

# ---- report ---------------------------------------------------------------------------------
# The MODE is recorded on the artefact itself: a green `fast` run must never be mistaken for a
# green full run, and an artefact that does not say which it was invites exactly that.
if (!NO_REPORT) {
  dir.create("output/report", recursive = TRUE, showWarnings = FALSE)
  mode_lab <- if (FAST) "fast (long suites skipped)" else "full"
  png_path <- "output/report/test_suite.png"
  ok_png <- FALSE
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    suppressMessages(library(ggplot2))
    d <- res; d$file <- factor(d$file, levels = rev(d$file))
    d$status <- factor(d$status, levels = c("PASS", "FAIL", "SKIP"))
    d$lab <- ifelse(d$status == "SKIP", "skipped",
                    ifelse(d$npass + d$nfail > 0, sprintf("%d/%d", d$npass, d$npass + d$nfail), d$note))
    p <- ggplot(d, aes(x = secs, y = file, fill = status)) +
      geom_col(width = 0.68) +
      geom_text(aes(label = lab), hjust = -0.15, size = 3.1, colour = "#333333") +
      scale_fill_manual(values = c(PASS = "#3f8f5b", FAIL = "#c0392b", SKIP = "#9aa0a6"), drop = FALSE) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
      labs(title = sprintf("gamble-core test suite - %s", if (nf == 0L) "ALL PASS" else sprintf("%d FILE(S) FAILED", nf)),
           subtitle = sprintf("%d files | %d assertions | %.1f min | mode: %s | %s",
                              nrow(res), sum(res$npass) + sum(res$nfail), sum(res$secs) / 60,
                              mode_lab, format(Sys.time(), "%Y-%m-%d %H:%M")),
           x = "seconds", y = NULL, fill = NULL) +
      theme_minimal(base_size = 11) +
      theme(panel.grid.major.y = element_blank(),
            plot.title = element_text(face = "bold"),
            legend.position = "top")
    ok_png <- tryCatch({
      ggsave(png_path, p, width = 9, height = max(3, 0.42 * nrow(res) + 1.6), dpi = 120)
      TRUE }, error = function(e) { cat("  [warn] could not write PNG: ", conditionMessage(e), "\n"); FALSE })
  } else cat("  [warn] ggplot2 not available -- HTML written without the figure\n")

  esc <- function(s) { s <- gsub("&", "&amp;", s); s <- gsub("<", "&lt;", s); gsub(">", "&gt;", s) }
  b64 <- function(p) if (file.exists(p)) {
    r <- readBin(p, "raw", file.info(p)$size)
    sprintf('<img src="data:image/png;base64,%s" alt="test suite results">', jsonlite::base64_enc(r))
  } else ""
  badge <- function(s) sprintf('<span class="b %s">%s</span>', tolower(s), s)
  rows <- paste(vapply(seq_len(nrow(res)), function(i) sprintf(
    '<tr class="%s"><td><code>%s</code></td><td>%s</td><td class="n">%.1fs</td><td class="n">%s</td><td>%s</td></tr>%s',
    tolower(res$status[i]), esc(res$file[i]), badge(res$status[i]), res$secs[i],
    if (res$npass[i] + res$nfail[i] > 0) sprintf("%d/%d", res$npass[i], res$npass[i] + res$nfail[i]) else "-",
    esc(res$note[i]),
    if (!is.null(fail_lines[[res$file[i]]]) && length(fail_lines[[res$file[i]]]))
      sprintf('<tr class="det"><td colspan="5"><pre>%s</pre></td></tr>',
              esc(paste(fail_lines[[res$file[i]]], collapse = "\n"))) else ""), character(1)), collapse = "\n")

  CSS <- "
:root{--ink:#1c1e21;--mut:#6b7280;--line:#e5e7eb;--pass:#3f8f5b;--fail:#c0392b;--skip:#9aa0a6}
*{box-sizing:border-box}body{margin:0;background:#f7f7f8;color:var(--ink);
 font:14px/1.5 -apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif}
header.top{background:#fff;border-bottom:1px solid var(--line);padding:22px 28px}
header.top h1{margin:0;font-size:19px}header.top p{margin:4px 0 0;color:var(--mut);font-size:13px}
.wrap{max-width:1040px;margin:0 auto;padding:22px 28px 60px}
.tiles{display:flex;gap:12px;flex-wrap:wrap;margin:0 0 20px}
.tile{background:#fff;border:1px solid var(--line);border-radius:8px;padding:12px 16px;min-width:120px}
.tval{font-size:22px;font-weight:600}.tlabel{color:var(--mut);font-size:12px;margin-top:2px}
figure{margin:0 0 22px;background:#fff;border:1px solid var(--line);border-radius:8px;padding:12px}
figure img{max-width:100%;display:block}
table{width:100%;border-collapse:collapse;background:#fff;border:1px solid var(--line);border-radius:8px;overflow:hidden}
th,td{padding:9px 12px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{background:#fafafa;font-size:12px;color:var(--mut);text-transform:uppercase;letter-spacing:.04em}
td.n{text-align:right;font-variant-numeric:tabular-nums}
tr.det td{background:#fff6f6;border-top:0}
pre{margin:0;font:12px/1.45 ui-monospace,SFMono-Regular,Menlo,monospace;white-space:pre-wrap;color:#7a1f16}
.b{display:inline-block;padding:2px 8px;border-radius:999px;font-size:11px;font-weight:600;color:#fff}
.b.pass{background:var(--pass)}.b.fail{background:var(--fail)}.b.skip{background:var(--skip)}
footer{color:var(--mut);font-size:12px;margin-top:18px}
code{font:12px ui-monospace,SFMono-Regular,Menlo,monospace}"

  html <- sprintf('<!doctype html><html><head><meta charset="utf-8">
<title>gamble-core test suite</title><style>%s</style></head><body>
<header class="top"><h1>gamble-core - test suite</h1>
<p>%s &middot; mode: <strong>%s</strong> &middot; generated %s</p></header>
<div class="wrap">
<div class="tiles">
 <div class="tile"><div class="tval">%d</div><div class="tlabel">files</div></div>
 <div class="tile"><div class="tval" style="color:var(--pass)">%d</div><div class="tlabel">passed</div></div>
 <div class="tile"><div class="tval" style="color:var(--fail)">%d</div><div class="tlabel">failed</div></div>
 <div class="tile"><div class="tval" style="color:var(--skip)">%d</div><div class="tlabel">skipped</div></div>
 <div class="tile"><div class="tval">%d</div><div class="tlabel">assertions</div></div>
 <div class="tile"><div class="tval">%.1f</div><div class="tlabel">minutes</div></div>
</div>
<figure>%s</figure>
<table><thead><tr><th>file</th><th>status</th><th>time</th><th>assertions</th><th>result</th></tr></thead>
<tbody>%s</tbody></table>
<footer>A file asserting nothing counts as a FAIL, not a pass. SKIP means the file declared its
inputs absent (<code>projects/</code> and <code>output/</code> are git-ignored, so project-dependent
tests only run where that build exists).%s</footer>
</div></body></html>',
    CSS, if (nf == 0L) "ALL PASS" else sprintf("%d FILE(S) FAILED", nf), mode_lab,
    format(Sys.time(), "%Y-%m-%d %H:%M"),
    nrow(res), sum(res$status == "PASS"), nf, ns, sum(res$npass) + sum(res$nfail), sum(res$secs) / 60,
    if (ok_png) b64(png_path) else "<p style='color:#6b7280;margin:8px'>figure unavailable</p>",
    rows,
    if (FAST) " <strong>This was a FAST run - the two long suites did not execute.</strong>" else "")
  writeLines(html, "output/report/test_suite.html")
  cat(sprintf("\n report: output/report/test_suite.html%s\n",
              if (ok_png) sprintf("  |  figure: %s", png_path) else ""))
}

if (nf > 0L) quit(status = 1L)
