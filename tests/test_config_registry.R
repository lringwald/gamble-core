#!/usr/bin/env Rscript
# =============================================================================
# test_config_registry.R — config/knobs.json is the one declaration; guard the derivations
# =============================================================================
# The generated artifacts (the shell mirror entrypoint.sh reads, the routine schema the platform
# renders) are only worth generating if they cannot silently drift from the registry. That drift is
# not hypothetical here: the hand-maintained routine enum left GLOBIOM_subclass unrunnable for
# weeks -- the code path existed, the enum entry did not, and nothing could notice.
# No JSON parser is used, deliberately: adding an R dependency to check a config file would be a
# worse trade than reading the two files as text.
# =============================================================================
np <- 0L; nf <- 0L
ok <- function(c_, m) { if (isTRUE(c_)) { np <<- np+1L; cat(sprintf("[PASS] %s\n", m)) }
                        else { nf <<- nf+1L; cat(sprintf("[FAIL] %s\n", m)) } }

for (f in c("config/knobs.json", "config/knobs.generated.sh",
            "docs/routine_config.schema.json", "docs/routine_config.schema.full.json")) {
  if (!file.exists(f)) { cat(sprintf("SKIP: %s missing -- run tools/gen_config.py\n", f)); quit(status = 0) }
}
reg_txt <- paste(readLines("config/knobs.json", warn = FALSE), collapse = "\n")
gen_txt <- paste(readLines("config/knobs.generated.sh", warn = FALSE), collapse = "\n")
sch_txt  <- paste(readLines("docs/routine_config.schema.json", warn = FALSE), collapse = "\n")
# Knob COVERAGE is asserted against the FULL schema. The lean one omits most knobs on purpose --
# a field declared there would hold a default that outranks the profile it should inherit from --
# so asserting coverage against it would only push the lean form back to twenty-five fields.
full_txt <- paste(readLines("docs/routine_config.schema.full.json", warn = FALSE), collapse = "\n")

reg_names <- regmatches(reg_txt, gregexpr('"name"\\s*:\\s*"[^"]+"', reg_txt))[[1]]
reg_names <- sub('.*"([^"]+)"$', "\\1", reg_names)
ok(length(reg_names) > 0, sprintf("registry declares %d knobs", length(reg_names)))

gen_names <- regmatches(gen_txt, gregexpr("KNOB_DEFAULT_[A-Z_0-9]+", gen_txt))[[1]]
gen_names <- sub("KNOB_DEFAULT_", "", gen_names)
ok(setequal(reg_names, gen_names),
   sprintf("shell mirror matches the registry (%d vs %d; missing: %s)", length(reg_names), length(gen_names),
           paste(c(setdiff(reg_names, gen_names), setdiff(gen_names, reg_names)), collapse = ", ")))

# The routine renders the form from a "root" key. A schema written at top level is not read at
# all, and the failure is silent -- an empty form, not an error.
ok(grepl('"root"\\s*:', sch_txt), "the routine schema is wrapped in its \"root\" envelope")
ok(grepl('"properties"\\s*:', sch_txt), "root carries a properties object")
ok(grepl('"root"\\s*:', full_txt), "the FULL schema carries the envelope too")
# The lean form must stay lean, and must not declare a knob whose default would beat the profile.
lean_fields <- regmatches(sch_txt, gregexpr('"[A-Z_]+"\\s*:\\s*\\{', sch_txt))[[1]]
lean_fields <- unique(sub('"([A-Z_]+)".*', "\\1", lean_fields))
ok(length(setdiff(lean_fields, c("TASK","PROFILE","CLASSIFICATION","VARIANT","RE_BLOCK","USE_BART","SUBSAMPLE","NITER","N_CHAINS","RUN_ID","GAMBLE_WORK_DIR","EXTRA"))) == 0,
   sprintf("the lean schema declares only the lean fields (%d)", length(lean_fields)))

miss_schema <- reg_names[!vapply(reg_names, function(n) grepl(sprintf('"%s"\\s*:', n), full_txt), TRUE)]
ok(length(miss_schema) == 0,
   sprintf("every registry knob reaches the FULL routine schema%s",
           if (length(miss_schema)) paste0(" -- absent: ", paste(miss_schema, collapse = ", ")) else ""))

# every key any profile sets must be declared, or the profile silently does nothing
pf <- list.files("config/profiles", pattern = "\\.env$", full.names = TRUE)
ok(length(pf) > 0, sprintf("%d profile(s) present", length(pf)))
bad <- character(0)
for (f in pf) {
  ln <- readLines(f, warn = FALSE)
  ln <- ln[nzchar(trimws(ln)) & !grepl("^\\s*#", ln)]
  k  <- trimws(sub("=.*$", "", ln))
  bad <- c(bad, sprintf("%s: %s", basename(f), setdiff(k, reg_names)[!is.na(setdiff(k, reg_names))]))
}
bad <- bad[!grepl(": $", bad)]
ok(length(bad) == 0,
   sprintf("no profile sets an undeclared key%s",
           if (length(bad)) paste0(" -- ", paste(bad, collapse = "; ")) else ""))

# the scope column is what makes a deviation announceable; it must be present and valid
scopes <- regmatches(reg_txt, gregexpr('"scope"\\s*:\\s*"[^"]+"', reg_txt))[[1]]
scopes <- sub('.*"([^"]+)"$', "\\1", scopes)
ok(length(scopes) == length(reg_names), "every knob declares a scope")
ok(all(scopes %in% c("model", "operational")), "every scope is model or operational")

cat(sprintf("RESULT: %d/%d checks passed\n", np, np + nf))
if (nf > 0L) quit(status = 1)
