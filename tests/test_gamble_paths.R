#!/usr/bin/env Rscript
# =============================================================================
# test_gamble_paths.R — the generated-artifact path resolver
# =============================================================================
# Guards the rule the layout depends on: output/ is regenerable and results/ holds fitted runs,
# both env-overridable so a container can place them on a mounted drive. The legacy lookup is the
# part most worth pinning -- readers must keep finding files written before the move, or a
# reorganisation silently orphans every run already on disk.
# =============================================================================
suppressMessages(source("codes/mnl_aux_func.R"))
np <- 0L; nf <- 0L
ok <- function(c_, m) { if (isTRUE(c_)) { np <<- np+1L; cat(sprintf("[PASS] %s\n", m)) }
                        else { nf <<- nf+1L; cat(sprintf("[FAIL] %s\n", m)) } }

.old <- Sys.getenv(c("GAMBLE_OUTPUT_DIR", "GAMBLE_RESULTS_DIR"), unset = NA)
T <- tempfile(); dir.create(T)
on.exit({ unlink(T, recursive = TRUE)
          for (n in names(.old)) if (is.na(.old[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, as.list(.old[n])) },
        add = TRUE)
Sys.setenv(GAMBLE_OUTPUT_DIR = file.path(T, "out"), GAMBLE_RESULTS_DIR = file.path(T, "res"))

ok(gamble_path("designs", "d.rds", create = FALSE) == file.path(T, "out", "designs", "d.rds"),
   "designs resolves under the output root")
ok(gamble_path("scratch", "g.rds", create = FALSE) == file.path(T, "out", "scratch", "g.rds"),
   "scratch resolves under the output root")
ok(gamble_run_path("PROJ", "prior_x", "lineage.txt", create = FALSE) ==
     file.path(T, "res", "gamble_model", "PROJ", "prior_x", "lineage.txt"),
   "a run resolves under results/gamble_model/<PROJECT>/ (the depth model_report.R reads)")

invisible(gamble_path("designs", "made.rds"))
ok(dir.exists(file.path(T, "out", "designs")), "the directory is created on demand, the file is not")
ok(!file.exists(file.path(T, "out", "designs", "made.rds")), "gamble_path does not create the file itself")
ok(inherits(tryCatch(gamble_path("nope", "x"), error = function(e) e), "error"),
   "an unknown kind is an error, not a silently wrong path")

# LEGACY LOOKUP: a file at the pre-move flat location must still be found.
invisible(file.create(file.path(T, "out", "dat_pixel_FULL_OLD.rds"))); Sys.sleep(1.1)
invisible(file.create(gamble_path("intermediate", "dat_pixel_FULL_NEW.rds")))
f <- gamble_find("intermediate", "dat_pixel_FULL_*.rds")
ok(length(f) == 2, "gamble_find sees BOTH the new and the legacy location")
ok(basename(f[1]) == "dat_pixel_FULL_NEW.rds", "newest file comes first")
ok(any(grepl("OLD", f)), "a file written before the move is not orphaned")
ok(length(gamble_find("designs", "nothing_matches_*.rds")) == 0, "no matches returns empty, not an error")

Sys.unsetenv(c("GAMBLE_OUTPUT_DIR", "GAMBLE_RESULTS_DIR"))
ok(gamble_output_root() == "output" && gamble_results_root() == "results",
   "defaults are output/ and results/ when the env is unset")

cat(sprintf("RESULT: %d/%d checks passed\n", np, np + nf))
if (nf > 0L) quit(status = 1)
