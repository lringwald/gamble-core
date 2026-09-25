# tools/gamble_config.R
# Pure R configuration compiler and validator for gamble-core container dispatch.

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Create a validated gamble-core job configuration
#'
#' @param task Character. Primary routine to run.
#' @param profile Character. Named profile from config/profiles/ or "none".
#' @param classification Character. Target variable to predict.
#' @param architecture Character or NULL. Inferred automatically from task if NULL.
#' @param variant Character. "factorized" or "iv" (nested only).
#' @param re_block Character. "intercept" or "intercept+socio" (nested only).
#' @param use_bart Logical. Whether to use BART for terrain (flat only).
#' @param design_path Character. Path to design dump or "auto".
#' @param subsample Integer or "auto". Subsample size.
#' @param nsample Integer or "auto". Kept draws per chain.
#' @param nburn Integer or "auto". Burn-in sweeps.
#' @param iterations Integer or "auto". Total sweeps (overrides nsample+nburn if set).
#' @param chains Integer or "auto". Number of parallel chains.
#' @param cores Integer or "auto". Number of worker cores.
#' @param run_id Character. Unique run tag or "auto".
#' @param work_dir Character. Mount path for outputs.
#' @param overrides Named list of key-value pairs for ad-hoc flags (checked against knobs.json).
#'
#' @return A validated list of configuration options.
gamble_config <- function(task = c("flat_fit", "nested", "flat_design", "report", 
                                   "recover_bart", "bart_gate", "count", "test"),
                          profile = "none",
                          classification = "GLOBIOM_subclass",
                          architecture = NULL,
                          variant = c("factorized", "iv"),
                          re_block = c("intercept", "intercept+socio"),
                          use_bart = FALSE,
                          design_path = "auto",
                          subsample = "auto",
                          nsample = "auto",
                          nburn = "auto",
                          iterations = "auto",
                          chains = "auto",
                          cores = "auto",
                          progress_sec = 600,
                          run_id = "auto",
                          work_dir = "/mnt/wdrv/gamble-core",
                          overrides = list()) {

  task <- match.arg(task)

  # 1. Infer architecture from task if omitted
  if (is.null(architecture)) {
    architecture <- if (task == "nested") "nested" else if (task == "count") "count" else "flat"
  } else {
    architecture <- match.arg(architecture, c("nested", "flat", "count", "utility"))
  }

  # 2. Cross-parameter validation
  if (architecture == "nested") {
    variant <- match.arg(variant)
    re_block <- match.arg(re_block)
    if (isTRUE(use_bart)) {
      stop("Configuration Error: 'use_bart = TRUE' is invalid for nested architecture. BART is flat-path only.")
    }
  }

  if (architecture == "flat") {
    if (!missing(variant) || !missing(re_block)) {
      warning("Configuration Notice: 'variant' and 're_block' are ignored when architecture is 'flat'.")
    }
  }

  # 3. Validate profile existence if not 'none'
  if (profile != "none" && file.exists("config/profiles")) {
    avail_profiles <- sub("\\.env$", "", list.files("config/profiles", pattern = "\\.env$"))
    if (!profile %in% avail_profiles) {
      warning(sprintf("Profile '%s' not found in config/profiles/. Available: %s",
                      profile, paste(avail_profiles, collapse = ", ")))
    }
  }

  # 4. Validate overrides against config/knobs.json if present
  if (length(overrides) > 0 && file.exists("config/knobs.json") && requireNamespace("jsonlite", quietly = TRUE)) {
    knobs_data <- jsonlite::fromJSON("config/knobs.json")
    declared_knobs <- knobs_data$knobs$name
    bad_keys <- setdiff(names(overrides), declared_knobs)
    if (length(bad_keys) > 0) {
      stop(sprintf("Configuration Error: Overrides contain undeclared knobs: %s. Container will reject them.",
                   paste(bad_keys, collapse = ", ")))
    }
  }

  list(
    execution = list(
      task = task,
      profile = profile,
      run_id = as.character(run_id)
    ),
    model = list(
      classification = classification,
      architecture = architecture,
      design_path = design_path,
      nested = if (architecture == "nested") list(variant = variant, re_block = re_block) else NULL,
      flat = if (architecture == "flat") list(use_bart = isTRUE(use_bart)) else NULL
    ),
    sampling = list(
      subsample = if (identical(subsample, "auto")) "auto" else as.integer(subsample),
      nsample = if (identical(nsample, "auto")) "auto" else as.integer(nsample),
      nburn = if (identical(nburn, "auto")) "auto" else as.integer(nburn),
      iterations = if (identical(iterations, "auto")) "auto" else as.integer(iterations),
      chains = if (identical(chains, "auto")) "auto" else as.integer(chains),
      cores = if (identical(cores, "auto")) "auto" else as.integer(cores),
      progress_sec = if (identical(progress_sec, "auto")) "auto" else as.integer(progress_sec)
    ),
    system = list(
      work_dir = work_dir
    ),
    overrides = overrides
  )
}

#' Flatten structured config into container environment variables
#'
#' @param cfg List produced by `gamble_config()` or loaded from JSON.
#' @return Named character vector ready for entrypoint.sh.
compile_to_env <- function(cfg) {
  env <- list()

  # 1. Execution
  env[["TASK"]] <- cfg$execution$task
  env[["PROFILE"]] <- cfg$execution$profile %||% "none"
  env[["RUN_ID"]] <- as.character(cfg$execution$run_id %||% "auto")

  # 2. Model Specification
  env[["CLASSIFICATION"]] <- cfg$model$classification %||% "GLOBIOM_subclass"
  if (!is.null(cfg$model$design_path) && cfg$model$design_path != "auto") {
    env[["DESIGN_PATH"]] <- cfg$model$design_path
  }

  arch <- cfg$model$architecture %||% "flat"
  if (identical(arch, "nested")) {
    nested_opts <- cfg$model$nested %||% list()
    env[["VARIANT"]] <- nested_opts$variant %||% "factorized"
    env[["RE_BLOCK"]] <- nested_opts$re_block %||% "intercept"
    env[["USE_BART"]] <- "FALSE"
  } else {
    flat_opts <- cfg$model$flat %||% list()
    env[["VARIANT"]] <- "factorized"
    env[["RE_BLOCK"]] <- "intercept"
    env[["USE_BART"]] <- if (isTRUE(flat_opts$use_bart)) "TRUE" else "FALSE"
  }

  # 3. Sampling
  sampling <- cfg$sampling %||% list()
  env[["SUBSAMPLE"]] <- as.character(sampling$subsample %||% "auto")
  if (!identical(sampling$nsample, "auto") && !is.null(sampling$nsample)) {
    env[["NSAMPLE"]] <- as.character(sampling$nsample)
  }
  if (!identical(sampling$nburn, "auto") && !is.null(sampling$nburn)) {
    env[["NBURN"]] <- as.character(sampling$nburn)
  }
  env[["NITER"]] <- as.character(sampling$iterations %||% "auto")
  env[["N_CHAINS"]] <- as.character(sampling$chains %||% "auto")
  env[["N_CORES"]] <- as.character(sampling$cores %||% "auto")
  if (!identical(sampling$progress_sec, "auto") && !is.null(sampling$progress_sec)) {
    env[["PROGRESS_SEC"]] <- as.character(sampling$progress_sec)
  }

  # 4. System Mounts
  system_opts <- cfg$system %||% list()
  env[["GAMBLE_WORK_DIR"]] <- system_opts$work_dir %||% "/mnt/wdrv/gamble-core"

  # 5. Overrides
  overrides <- cfg$overrides %||% list()
  if (length(overrides) > 0) {
    extra_pairs <- vapply(names(overrides), function(k) {
      paste0(k, "=", overrides[[k]])
    }, FUN.VALUE = character(1))
    env[["EXTRA"]] <- paste(extra_pairs, collapse = ", ")
  } else {
    env[["EXTRA"]] <- "none"
  }

  stats::setNames(as.character(unlist(env)), names(env))
}

#' Export environment variables to file, container flags, or JSON for Accelerator UI
#'
#' @param env_vars Named character vector from `compile_to_env()`.
#' @param format Target format: "dotenv", "bash", "docker", "singularity", or "json".
#' @param filepath Optional destination path to write output.
#' @return Formatted character string or lines (invisible if written to file).
export_env <- function(env_vars, 
                       format = c("dotenv", "bash", "docker", "singularity", "json"), 
                       filepath = NULL) {
  format <- match.arg(format)

  lines <- switch(format,
    "dotenv" = paste0(names(env_vars), "=", env_vars),
    "bash" = paste0("export ", names(env_vars), "=\"", env_vars, "\""),
    "docker" = paste(paste0("-e ", names(env_vars), "=\"", env_vars, "\""), collapse = " "),
    "singularity" = paste(paste0("--env ", names(env_vars), "=\"", env_vars, "\""), collapse = " "),
    "json" = {
      if (requireNamespace("jsonlite", quietly = TRUE)) {
        jsonlite::toJSON(as.list(env_vars), pretty = TRUE, auto_unbox = TRUE)
      } else {
        paste0("{\n", paste(sprintf('  "%s": "%s"', names(env_vars), env_vars), collapse = ",\n"), "\n}")
      }
    }
  )

  if (!is.null(filepath)) {
    writeLines(lines, con = filepath)
    message(sprintf("Configuration exported successfully to '%s' (%s format).", filepath, format))
  }

  invisible(lines)
}

#' Load and compile a JSON job configuration
#'
#' @param json_path Path to JSON file.
#' @return Named character vector of environment variables.
load_and_compile_json <- function(json_path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to parse JSON config files.")
  }

  raw_cfg <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)
  
  # Validate model-specific exclusivity
  arch <- raw_cfg$model$architecture
  if (identical(arch, "nested") && isTRUE(raw_cfg$model$flat$use_bart)) {
    stop("Validation failed: 'use_bart' cannot be set to TRUE on a nested architecture.")
  }

  compile_to_env(raw_cfg)
}
