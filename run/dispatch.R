# run/dispatch.R
# Interactive workflow script for compiling and dispatching gamble-core configurations.
# Source tools/gamble_config.R to configure, validate, and export container jobs.

source("tools/gamble_config.R")

cat("======================================================================\n")
cat(" gamble-core Routine Dispatch Compiler\n")
cat("======================================================================\n\n")

# -------------------------------------------------------------------------
# Scenario A: Fast Smoke Test (Flat Fit on 5k Subsample)
# -------------------------------------------------------------------------
cat(">>> Scenario A: Fast Smoke Test (Flat Fit)\n")
smoke_job <- gamble_config(
  task = "flat_fit",
  classification = "GLOBIOM_subclass",
  subsample = 5000,
  nsample = 50,
  nburn = 20,
  chains = 1,
  cores = 1,
  run_id = "smoke_test_01"
)

smoke_env <- compile_to_env(smoke_job)
print(smoke_env)

# Export to a local .env file
export_env(smoke_env, format = "dotenv", filepath = ".env.smoke")
cat("\n")

# -------------------------------------------------------------------------
# Scenario B: Full Production Nested Model (Factorized Nests)
# -------------------------------------------------------------------------
cat(">>> Scenario B: Full Production Nested Model\n")
nested_job <- gamble_config(
  task = "nested",
  profile = "globiom_nested_prod",
  classification = "GLOBIOM_subclass",
  variant = "factorized",
  re_block = "intercept",
  chains = 4,
  cores = 4,
  run_id = "prod_nested_2026"
)

nested_env <- compile_to_env(nested_job)

# Export as JSON ready to paste into the Accelerator Routine UI JSON editor:
cat("--- Accelerator JSON Editor Payload ---\n")
cat(export_env(nested_env, format = "json"), "\n\n")

# -------------------------------------------------------------------------
# Scenario C: Flat Production Fit with BART on Terrain
# -------------------------------------------------------------------------
cat(">>> Scenario C: Flat Model with BART (Command Generation)\n")
prod_bart_job <- gamble_config(
  task = "flat_fit",
  profile = "globiom_flat_bart",
  classification = "GLOBIOM_subclass",
  use_bart = TRUE,
  nsample = 6000,
  nburn = 2000,
  chains = 4,
  cores = 4,
  run_id = "prod_bart_run",
  overrides = list(
    BART_COLS = "topo",
    NTREES_BART = "25"
  )
)

prod_bart_env <- compile_to_env(prod_bart_job)

# Generate Docker CLI flags
docker_flags <- export_env(prod_bart_env, format = "docker")
cat("--- Docker CLI Command ---\n")
cat("docker run --rm", docker_flags, "gamble-core:latest\n\n")

# Generate Apptainer / Singularity CLI flags
sing_flags <- export_env(prod_bart_env, format = "singularity")
cat("--- Apptainer / Singularity CLI Command ---\n")
cat("apptainer run", sing_flags, "gamble-core.sif\n\n")

# -------------------------------------------------------------------------
# Scenario D: Direct Local Execution Wrapper
# -------------------------------------------------------------------------
dispatch_locally <- function(env_vars, entrypoint_path = "./entrypoint.sh") {
  if (!file.exists(entrypoint_path)) {
    stop("entrypoint script not found at ", entrypoint_path)
  }
  do.call(Sys.setenv, as.list(env_vars))
  message("Environment exported. Executing ", entrypoint_path, "...")
  status <- system(entrypoint_path)
  if (status != 0) {
    stop(sprintf("Task execution failed with exit code %d", status))
  }
  message("Task completed successfully.")
}

cat("To run locally with the active R session environment:\n")
cat("  dispatch_locally(smoke_env)\n")
