# tests/test_gamble_config.R
# Test suite for gamble_config compiler and validator

source("tools/gamble_config.R")

test_that("gamble_config infers architecture correctly", {
  cfg_flat <- gamble_config(task = "flat_fit")
  expect_equal(cfg_flat$model$architecture, "flat")

  cfg_nested <- gamble_config(task = "nested")
  expect_equal(cfg_nested$model$architecture, "nested")

  cfg_count <- gamble_config(task = "count")
  expect_equal(cfg_count$model$architecture, "count")
})

test_that("gamble_config blocks invalid cross-parameter combinations", {
  expect_error(
    gamble_config(task = "nested", use_bart = TRUE),
    "use_bart = TRUE.*invalid for nested architecture"
  )
})

test_that("gamble_config rejects undeclared overrides", {
  expect_error(
    gamble_config(task = "flat_fit", overrides = list(TOTALLY_FAKE_KNOB_XYZ = "123")),
    "Overrides contain undeclared knobs"
  )
})

test_that("compile_to_env flattens variables properly", {
  cfg <- gamble_config(
    task = "flat_fit",
    classification = "GLOBIOM_subclass",
    nsample = 6000,
    nburn = 2000,
    chains = 4,
    cores = 4,
    use_bart = TRUE,
    overrides = list(BART_COLS = "topo")
  )
  env <- compile_to_env(cfg)

  expect_equal(env[["TASK"]], "flat_fit")
  expect_equal(env[["USE_BART"]], "TRUE")
  expect_equal(env[["NSAMPLE"]], "6000")
  expect_equal(env[["NBURN"]], "2000")
  expect_equal(env[["N_CHAINS"]], "4")
  expect_equal(env[["N_CORES"]], "4")
  expect_equal(env[["EXTRA"]], "BART_COLS=topo")
})

test_that("export_env produces valid formats", {
  cfg <- gamble_config(task = "flat_fit", run_id = "test_01")
  env <- compile_to_env(cfg)

  json_out <- export_env(env, format = "json")
  expect_true(grepl('"TASK": "flat_fit"', json_out))

  docker_out <- export_env(env, format = "docker")
  expect_true(grepl('-e TASK="flat_fit"', docker_out))

  sing_out <- export_env(env, format = "singularity")
  expect_true(grepl('--env TASK="flat_fit"', sing_out))

  bash_out <- export_env(env, format = "bash")
  expect_true(any(grepl('^export TASK="flat_fit"', bash_out)))
})

cat("All gamble_config tests passed successfully!\n")
