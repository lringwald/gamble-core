# Consolidate the BOV + SGT composition δ-fits into one tidy parameter table:
# species × driver × subclass (D/O/F) -> δ mean + 90% CI + visible flag.
suppressMessages(library(data.table))
rows <- list()
for (sp in c("bov","sgt")) {
  f <- file.path("output/composition", paste0(sp, "_composition_fit.rds"))
  if (!file.exists(f)) next
  fit <- readRDS(f); dr <- fit$drivers
  for (cl in c("D","O","F")) {
    M <- fit[[paste0("d", cl)]]                         # [k x draws] sum-to-zero δ
    mean_ <- rowMeans(M); lo <- apply(M,1,quantile,.05); hi <- apply(M,1,quantile,.95)
    rows[[paste(sp,cl)]] <- data.table(species=toupper(sp), subclass=paste0(toupper(sp),cl),
      driver=dr, delta=mean_, lo=lo, hi=hi, visible=(sign(lo)==sign(hi)))
  }
}
tab <- rbindlist(rows)[driver != "intercept"]
fwrite(tab, "output/composition/subclass_allocation_parameters.csv")
cat(sprintf("Wrote %d rows (6 subclasses) -> output/composition/subclass_allocation_parameters.csv\n", nrow(tab)))
cat("\nVisible drivers per subclass:\n")
print(tab[, .(n_drivers=.N, n_visible=sum(visible)), by=subclass])
cat("\nStrongest visible driver per subclass:\n")
print(tab[visible == TRUE][order(subclass, -abs(delta))][, .SD[1], by=subclass][, .(subclass, driver, delta=round(delta,2))])
