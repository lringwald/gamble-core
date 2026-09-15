# Consolidate the BOV + SGT composition δ-fits into one tidy parameter table:
# species × driver × subclass -> δ mean + 90% CI + visible flag.
#
# CATEGORIES COME FROM THE FIT, not from a hardcoded c("D","O","F"): a species may be fitted on the
# 3-category D/O/F target or on a finer one (the 8-category CAPRI cattle set), and the two can differ
# BETWEEN species in the same run. A fixed list silently produced NULL deltas for anything that was
# not D/O/F -- the table still wrote, one row block short.
suppressMessages(library(data.table))
rows <- list()
for (sp in c("bov","sgt")) {
  f <- file.path("output/composition", paste0(sp, "_composition_fit.rds"))
  if (!file.exists(f)) next
  fit <- readRDS(f); dr <- fit$drivers
  # general form first; fall back to the legacy dD/dO/dF fields for fits saved before 2026-09-14
  cats <- if (!is.null(fit$cats)) fit$cats else c("D","O","F")
  getd <- function(cl) if (!is.null(fit$delta)) fit$delta[[cl]] else fit[[paste0("d", cl)]]
  cat(sprintf(">>> %s: %d categories (%s)\n", sp, length(cats), paste(cats, collapse = ", ")))
  for (cl in cats) {
    M <- getd(cl)                                       # [k x draws] sum-to-zero δ
    if (is.null(M)) { warning(sprintf("%s: no delta for category '%s' -- skipped", sp, cl), call. = FALSE); next }
    mean_ <- rowMeans(M); lo <- apply(M,1,quantile,.05); hi <- apply(M,1,quantile,.95)
    rows[[paste(sp,cl)]] <- data.table(species=toupper(sp), subclass=paste0(toupper(sp),cl),
      driver=dr, delta=mean_, lo=lo, hi=hi, visible=(sign(lo)==sign(hi)))
  }
}
if (!length(rows)) stop("no composition fits found under output/composition/ -- run fit_composition.R first")
tab <- rbindlist(rows)[driver != "intercept"]
fwrite(tab, "output/composition/subclass_allocation_parameters.csv")
cat(sprintf("Wrote %d rows (%d subclasses) -> output/composition/subclass_allocation_parameters.csv\n",
            nrow(tab), uniqueN(tab$subclass)))
cat("\nVisible drivers per subclass:\n")
print(tab[, .(n_drivers=.N, n_visible=sum(visible)), by=subclass])
cat("\nStrongest visible driver per subclass:\n")
print(tab[visible == TRUE][order(subclass, -abs(delta))][, .SD[1], by=subclass][, .(subclass, driver, delta=round(delta,2))])
