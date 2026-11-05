#!/usr/bin/env Rscript
# ==============================================================================
# scripts/generate_manifest.R
# Generates Atomic Task Manifests: config/manifest_synth.csv & config/manifest_real.csv
# Supports --smoke_test flag for local end-to-end validation
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
})

source("R/01_ode_simulators.R")
source("R/10_real_data_experiments.R")

args <- commandArgs(trailingOnly = TRUE)
is_smoke <- "--smoke_test" %in% args

cfg <- yaml::read_yaml("config/benchmark_grid.yaml")
methods <- cfg$methods
sys_reg <- get_systems_registry()
real_reg <- get_real_datasets_registry()

if (is_smoke) {
  message("[MANIFEST] Generating SMOKE TEST manifests (--smoke_test)...")
  synth_systems <- cfg$smoke_test$synthetic_systems
  synth_Ns <- cfg$smoke_test$synthetic_sample_sizes
  synth_noises <- cfg$smoke_test$synthetic_noise_eta
  synth_reps <- seq_len(cfg$smoke_test$n_replications)

  real_datasets <- cfg$smoke_test$real_datasets
  real_diffs <- cfg$smoke_test$real_diff_methods
  real_reps <- seq_len(cfg$smoke_test$n_replications)
  fast_flag <- 1L
} else {
  message("[MANIFEST] Generating FULL ICML benchmark manifests...")
  synth_systems <- cfg$synthetic_experiment$systems
  synth_Ns <- cfg$synthetic_experiment$sample_sizes
  synth_noises <- cfg$synthetic_experiment$derivative_noise_eta
  synth_reps <- seq_len(cfg$synthetic_experiment$n_replications)

  real_datasets <- cfg$real_data_experiment$datasets
  real_diffs <- cfg$real_data_experiment$diff_methods
  real_reps <- seq_len(cfg$real_data_experiment$n_replications)
  fast_flag <- 0L
}

# 1. Build Synthetic Manifest
synth_rows <- list()
for (s_idx in seq_along(synth_systems)) {
  sys_key <- synth_systems[s_idx]
  sys_num <- sys_reg[[sys_key]]$sys_id
  sys_grp <- sys_reg[[sys_key]]$group
  for (n_idx in seq_along(synth_Ns)) {
    n_val <- synth_Ns[n_idx]
    for (eta_idx in seq_along(synth_noises)) {
      eta_val <- synth_noises[eta_idx]
      for (rep_id in synth_reps) {
        # Deterministic seed formula from Section 2.2:
        # seed = 10000 * sys_id + 1000 * noise_id + 100 * n_id + rep_id
        seed_val <- as.integer(10000L * sys_num + 1000L * eta_idx + 100L * n_idx + rep_id)
        for (m_id in methods) {
          synth_rows[[length(synth_rows) + 1L]] <- list(
            suite = "synth",
            system_id = sys_key,
            sys_num = sys_num,
            system_group = sys_grp,
            sample_size = as.integer(n_val),
            noise_eta = as.numeric(eta_val),
            noise_regime = "derivative_snr",
            rep_id = as.integer(rep_id),
            seed = seed_val,
            method_id = m_id,
            fast_mode = fast_flag
          )
        }
      }
    }
  }
}
dt_synth <- rbindlist(synth_rows)
dt_synth[, task_id := .I]
setcolorder(dt_synth, c("task_id", setdiff(names(dt_synth), "task_id")))
fwrite(dt_synth, "config/manifest_synth.csv")
message("[OK] Wrote config/manifest_synth.csv with ", nrow(dt_synth), " atomic tasks.")

# 2. Build Real-World Manifest
real_rows <- list()
for (d_idx in seq_along(real_datasets)) {
  ds_key <- real_datasets[d_idx]
  ds_num <- real_reg[[ds_key]]$dataset_id
  ds_dom <- real_reg[[ds_key]]$domain
  for (dm_idx in seq_along(real_diffs)) {
    dm_val <- real_diffs[dm_idx]
    for (fold_id in real_reps) {
      seed_val <- as.integer(50000L + 1000L * ds_num + 100L * dm_idx + fold_id)
      for (m_id in methods) {
        real_rows[[length(real_rows) + 1L]] <- list(
          suite = "real",
          dataset_id = ds_key,
          dataset_num = ds_num,
          domain = ds_dom,
          diff_method = dm_val,
          fold_id = as.integer(fold_id),
          seed = seed_val,
          method_id = m_id,
          fast_mode = fast_flag
        )
      }
    }
  }
}
dt_real <- rbindlist(real_rows)
dt_real[, task_id := .I]
setcolorder(dt_real, c("task_id", setdiff(names(dt_real), "task_id")))
fwrite(dt_real, "config/manifest_real.csv")
message("[OK] Wrote config/manifest_real.csv with ", nrow(dt_real), " atomic tasks.")

