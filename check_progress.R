#!/usr/bin/env Rscript
# =============================================================================
# check_progress.R — how far along is a streamed sampler run?
# =============================================================================
#   Rscript check_progress.R                       # newest run under results/ or output/
#   Rscript check_progress.R <dir>                 # a run dir, or its posterior/ dir
#   Rscript check_progress.R <dir> --watch         # refresh every 60s until done
#
# Reads the batch files the sampler streams to disk. Note the two things that make a naive look at
# the directory misleading:
#   * NOTHING is written during burn-in, and then only every posterior_batch_size retained draws
#     (50 by default = 400 sweeps at thin 8). A quiet directory early on is normal, not a stall.
#   * With chains in future::multisession, the workers' console output does NOT reach the parent's
#     log, so the log going silent after the header says nothing about progress either.
# The reliable signals are the batch files' iteration stamps and their mtimes.
# =============================================================================
suppressMessages({library(qs2)})
args <- commandArgs(trailingOnly = TRUE)
WATCH <- "--watch" %in% args; args <- setdiff(args, "--watch")

find_latest <- function() {
  cand <- unlist(lapply(c("results/gamble_model", "output"), function(r)
    if (dir.exists(r)) list.dirs(r, recursive = TRUE) else character(0)))
  cand <- cand[basename(cand) == "posterior"]
  if (!length(cand)) return("")
  cand[order(file.mtime(cand), decreasing = TRUE)][1]
}
d <- if (length(args)) args[1] else find_latest()
if (!nzchar(d)) stop("no run found; pass a directory")
if (dir.exists(file.path(d, "posterior"))) d <- file.path(d, "posterior")
if (!dir.exists(d)) stop("not a directory: ", d)

report <- function(d) {
  meta <- tryCatch(qs_read(file.path(d, "model_metadata.qs")), error = function(e) NULL)
  lin  <- {
    lf <- file.path(dirname(d), "lineage.txt")
    if (file.exists(lf)) readLines(lf) else character(0)
  }
  # target sweeps: metadata first, then the lineage line written by estimate_prior.R
  niter <- meta$niter; nburn <- meta$nburn; thin <- meta$thin
  if (is.null(niter) && length(lin)) {
    m <- regmatches(lin, regexec("sweeps\\s+(\\d+) \\(burn (\\d+), thin (\\d+)\\)", lin))
    m <- m[lengths(m) == 4]
    if (length(m)) { niter <- as.integer(m[[1]][2]); nburn <- as.integer(m[[1]][3]); thin <- as.integer(m[[1]][4]) }
  }
  f <- list.files(d, "^posterior_batch_[0-9]+_chain_.*\\.qs$", full.names = TRUE)
  cat(sprintf("\n%s\n%s\n", d, strrep("-", min(nchar(d), 78))))
  if (length(lin)) cat(sprintf("  %s\n", paste(grep("^(segment|parent)", lin, value = TRUE), collapse = " | ")))

  # Heartbeats: the live signal, written every 100 sweeps INCLUDING during burn-in. Batch files only
  # appear later and much more rarely, so prefer this when present.
  hb <- list.files(d, "^progress_chain_.*\\.txt$", full.names = TRUE)
  if (length(hb)) {
    cat(sprintf("\n  %-7s %9s %9s %-12s %10s\n", "chain", "sweep", "of", "phase", "quiet"))
    it_all <- integer(0)
    for (h in sort(hb)) {
      ln <- tryCatch(strsplit(readLines(h, warn = FALSE)[1], " ")[[1]], error = function(e) NULL)
      if (is.null(ln) || length(ln) < 3) next
      it_all <- c(it_all, as.integer(ln[1]))
      cat(sprintf("  %-7s %9s %9s %-12s %8.1fm\n", sub(".*_chain_(.*)\\.txt$", "\\1", basename(h)),
          ln[1], ln[2], ln[3], as.numeric(difftime(Sys.time(), file.mtime(h), units = "mins"))))
    }
    if (length(it_all) && !is.null(niter)) {
      st <- if (!is.null(meta$started_at)) as.POSIXct(meta$started_at) else min(file.mtime(hb))
      el <- as.numeric(difftime(Sys.time(), st, units = "mins"))
      dn <- min(it_all)
      cat(sprintf("\n  slowest chain at sweep %d of %d (%.1f%%) after %.0f min\n", dn, niter, 100*dn/niter, el))
      if (dn > 0 && el > 0) cat(sprintf("  ~%.1f sweeps/min -> ETA %.1f h\n", dn/el, (niter-dn)/(dn/el)/60))
    }
  }
  if (!length(f)) {
    cat("  no batches yet.\n")
    if (!is.null(nburn)) cat(sprintf("  burn-in is %d sweeps and writes NOTHING; first batch lands after burn + %d sweeps.\n",
                                     nburn, (if (is.null(meta$posterior_batch_size)) 50L else meta$posterior_batch_size) * (thin %||% 1L)))
    return(invisible(NULL))
  }
  ch <- sub(".*_chain_(.*)\\.qs$", "\\1", basename(f))
  now <- Sys.time()
  rows <- lapply(sort(unique(ch)), function(cc) {
    fc <- f[ch == cc]; fc <- fc[order(as.integer(sub(".*batch_([0-9]+)_.*", "\\1", basename(fc))))]
    b  <- tryCatch(qs_read(fc[length(fc)]), error = function(e) NULL)   # last batch may be mid-write
    if (is.null(b)) { fc <- fc[-length(fc)]; if (!length(fc)) return(NULL); b <- qs_read(fc[length(fc)]) }
    it <- b[[length(b)]]$iter
    ll <- b[[length(b)]]$log_lik
    nd <- (length(fc) - 1L) * (if (is.null(meta$posterior_batch_size)) 50L else meta$posterior_batch_size) + length(b)
    data.frame(chain = cc, batches = length(fc), draws = nd, iter = if (is.null(it)) NA_integer_ else it,
               loglik = if (is.null(ll)) NA_real_ else ll,
               quiet_min = as.numeric(difftime(now, file.mtime(fc[length(fc)]), units = "mins")))
  })
  r <- do.call(rbind, Filter(Negate(is.null), rows))
  cat(sprintf("  %-7s %8s %7s %9s %12s %10s\n", "chain", "batches", "draws", "sweep", "log-lik", "quiet"))
  for (i in seq_len(nrow(r))) cat(sprintf("  %-7s %8d %7d %9s %12.0f %8.1fm\n",
      r$chain[i], r$batches[i], r$draws[i], if (is.na(r$iter[i])) "?" else format(r$iter[i]),
      r$loglik[i], r$quiet_min[i]))
  if (!is.null(niter) && !all(is.na(r$iter))) {
    done <- min(r$iter, na.rm = TRUE)
    # Rate over the WHOLE run, anchored on when it started. Measuring the span between batch mtimes
    # instead is worthless until several batches exist: with one batch per chain the span is seconds,
    # which reported ~991 sweeps/min and an ETA of 0.0 h for a run with two hours left.
    st <- if (!is.null(meta$started_at)) as.POSIXct(meta$started_at) else {
      bl <- grep("^built", lin, value = TRUE)
      if (length(bl)) as.POSIXct(trimws(sub("^built\\s+", "", bl[1]))) else min(file.mtime(f))
    }
    # `done` is a FLOOR that only advances when a batch lands, every posterior_batch_size*thin sweeps.
    # So measure the rate as of the batch that reported it, not as of now -- otherwise the elapsed
    # time keeps growing against a frozen sweep count and the ETA drifts upward between batches,
    # which reads as the run slowing down when nothing has changed.
    t_at <- max(file.mtime(f))
    el   <- as.numeric(difftime(t_at, st, units = "mins"))
    since<- as.numeric(difftime(Sys.time(), t_at, units = "mins"))
    per  <- if (el > 0) done / el else NA_real_
    step <- (if (is.null(meta$posterior_batch_size)) 50L else meta$posterior_batch_size) * (thin %||% 1L)
    cat(sprintf("\n  slowest chain at sweep >= %d of %d  (>= %.1f%%) as of %.0f min in\n",
                done, niter, 100 * done / niter, el))
    if (is.finite(per) && per > 0) {
      est <- min(niter, done + per * since)
      cat(sprintf("  ~%.1f sweeps/min -> now approx sweep %.0f (%.1f%%), ETA %.1f h\n",
                  per, est, 100 * est / niter, max(0, (niter - est) / per / 60)))
    }
    cat(sprintf("  the sweep counter advances only every %d sweeps (next update at %d) -- a frozen\n  number between batches is NOT a stall; judge liveness by CPU below\n",
                step, min(niter, done + step)))
    if (max(r$quiet_min) > 30)
      cat(sprintf("  NOTE: chain %s has not written for %.0f min -- check it is alive (see below)\n",
                  r$chain[which.max(r$quiet_min)], max(r$quiet_min)))
  }
  invisible(r)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

repeat {
  res <- report(d)
  cat("\n  live R processes (Rscript execs into R, so grep for the latter):\n")
  ps <- suppressWarnings(system("ps -A -o pid,etime,%cpu,command | grep '[R] --no-echo' | awk '{printf \"    pid %s  up %s  cpu %s%%\\n\", $1, $2, $3}'", intern = TRUE))
  cat(if (length(ps)) paste(ps, collapse = "\n") else "    none — the run is finished or dead", "\n\n")
  if (!WATCH) break
  Sys.sleep(60)
}
