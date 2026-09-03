#!/usr/bin/env Rscript
# =============================================================================
# run_all.R — one entry point for every test in the repo
# =============================================================================
# Why this exists: there was no single command that ran everything, so the pixel suite went
# unverified for a full day of sampler changes while targeted tests passed. Each file below is
# self-contained and exits non-zero on failure; this just runs them and prints one table.
#
#   Rscript tests/run_all.R            # everything
#   Rscript tests/run_all.R fast       # skip the two long suites
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
FAST <- "fast" %in% args
root <- if (basename(getwd()) == "tests") ".." else "."
setwd(root)

# NOTE: this list is hand-maintained, which is itself a hazard -- three feature tests
# (re_mean_shift, joint_shrink, mundlak) sat outside it and never ran here. When adding a
# test file, add it BELOW.
files <- c(
  "tests/test_nested_cut_units.R",    # fast helper units (M=1 IV fields) -- runs first, fails fast
  "tests/test_altspec_sampler.R",     # delta: recovery, multi-block, per_class, symmetric, scaling
  "tests/test_altspec_nested.R",      # blocks through nested_cut: wiring, IV path, interactions
  "tests/test_nested_cut.R",
  "tests/test_nested_cut_focal.R",
  "tests/test_nested_iv.R",
  "tests/test_re_mean_shift.R",       # sum-to-zero RE identification (mu == population average)
  "tests/test_joint_shrink.R",        # joint FE/RE gate kappa_v, incl. the seed-noise floor
  "tests/test_mundlak.R"              # group-mean prior mean: gamma recovery + beta debiasing
)
.known <- files
.found <- sort(list.files("tests", pattern = "^test_.*\\.R$", full.names = TRUE))
.miss  <- setdiff(.found, .known)
if (length(.miss))
  cat(sprintf("  [warn] test file(s) present but NOT in the run list: %s\n",
              paste(basename(.miss), collapse = ", ")))
if (!FAST) files <- c(files, "codes/test_suite_lu_pixel.R", "codes/test_suite_ls_count.R")
files <- files[file.exists(files)]

cat(sprintf("\n%s\n running %d test file(s)%s\n%s\n", strrep("=", 64), length(files),
            if (FAST) "  [fast: long suites skipped]" else "", strrep("=", 64)))
res <- data.frame(file = character(), status = character(), secs = numeric(), note = character())
for (f in files) {
  t0 <- Sys.time()
  out <- suppressWarnings(system2("Rscript", f, stdout = TRUE, stderr = TRUE))
  rc <- attr(out, "status"); rc <- if (is.null(rc)) 0L else as.integer(rc)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # a suite reports "RESULT: a/b"; a targeted test reports [PASS]/[FAIL] lines
  hit <- grep("RESULT:", out, value = TRUE)
  npass <- sum(grepl("\\[PASS\\]", out)); nfail <- sum(grepl("\\[FAIL\\]", out))
  # A file that asserts NOTHING is not a pass. Exit code 0 with zero assertions means it did not run
  # (an unreadable file, an early return, a stub) -- silent success is the failure mode this whole
  # repo keeps hitting, and a green runner that hides it is worse than no runner.
  silent <- !length(hit) && npass == 0L && nfail == 0L
  note <- if (length(hit)) trimws(sub(".*RESULT:", "", hit[1]))
          else if (silent) "NO ASSERTIONS RUN -- did the file execute?"
          else sprintf("%d pass / %d fail", npass, nfail)
  bad <- rc != 0L || nfail > 0L || silent || (length(hit) && grepl("FAILED", hit[1]))
  res <- rbind(res, data.frame(file = basename(f), status = if (bad) "FAIL" else "PASS",
                               secs = el, note = note))
  cat(sprintf("  [%s] %-34s %6.1fs  %s\n", if (bad) "FAIL" else "PASS", basename(f), el, note))
  if (bad) cat(paste0("        ", tail(grep("\\[FAIL\\]|Error", out, value = TRUE), 3), collapse = "\n"), "\n")
}
nf <- sum(res$status == "FAIL")
cat(sprintf("\n%s\n %s  (%d file(s), %.1f min)\n%s\n", strrep("=", 64),
            if (nf == 0L) "ALL PASS" else sprintf("%d FILE(S) FAILED", nf),
            nrow(res), sum(res$secs)/60, strrep("=", 64)))
if (nf > 0L) quit(status = 1L)
