# =============================================================================
# helper_recovery.R — capture truth-vs-estimate so a human can SEE recovery
# =============================================================================
# The suite's figure was a bar chart of runtime coloured by pass/fail: it showed which files ran,
# never whether anything was recovered. A [PASS] only says an assertion's threshold was met, and
# the reader has to trust that the threshold was chosen well.
#
# The tests already COMPUTE truth and estimate -- they print lines like "true delta: 0.800" and
# then discard them. tv_record() appends them to one CSV instead, and run_all.R renders the lot:
# every parameter against its truth, with its interval, grouped by feature. A reader can then check
# the claim rather than the count.
#
# Deliberately dependency-free and append-only: a test that crashes still leaves the rows it got to.
TV_OUT <- Sys.getenv("TESTVIZ_OUT", "output/report/test_recovery.csv")

tv_record <- function(feature, param, truth, est, lo = NA_real_, hi = NA_real_, note = "") {
  dir.create(dirname(TV_OUT), recursive = TRUE, showWarnings = FALSE)
  new <- !file.exists(TV_OUT)
  row <- data.frame(feature = feature, param = param, truth = as.numeric(truth),
                    est = as.numeric(est), lo = as.numeric(lo), hi = as.numeric(hi),
                    note = note, stringsAsFactors = FALSE)
  utils::write.table(row, TV_OUT, sep = ",", row.names = FALSE, col.names = new,
                     append = !new, qmethod = "double")
  invisible(row)
}

# Convenience for a whole vector of parameters (per-class effects, per-group REs, ...).
tv_record_vec <- function(feature, truth, est, lo = NULL, hi = NULL, names = NULL, note = "") {
  n <- length(truth); names <- if (is.null(names)) paste0("p", seq_len(n)) else names
  for (i in seq_len(n))
    tv_record(feature, names[i], truth[i], est[i],
              if (is.null(lo)) NA_real_ else lo[i], if (is.null(hi)) NA_real_ else hi[i], note)
  invisible(NULL)
}
