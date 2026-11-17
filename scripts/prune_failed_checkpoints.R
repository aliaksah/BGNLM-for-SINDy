#!/usr/bin/env Rscript
# ==============================================================================
# scripts/prune_failed_checkpoints.R
# Removes ONLY checkpoints with Status != "OK" (or corrupted/unreadable RDS files)
# Preserves all valid Status == "OK" checkpoints intact for seamless re-execution.
# ==============================================================================

if (file.exists("R/00_env_setup.R")) source("R/00_env_setup.R")

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

prune_suite <- function(suite = c("synth", "real"), dry_run = FALSE) {
  suite <- match.arg(suite)
  ckpt_dir <- paste0("results/checkpoints_", suite)
  files <- list.files(ckpt_dir, pattern = "^task_[0-9]+\\.rds$", full.names = TRUE)
  if (length(files) == 0L) {
    message("[PRUNE] No checkpoints found in ", ckpt_dir)
    return(invisible(NULL))
  }

  n_cores <- min(8L, max(1L, detectCores() - 1L))
  is_failed <- unlist(mclapply(files, function(f) {
    ok <- tryCatch({
      dt <- readRDS(f)
      !is.null(dt) && nrow(dt) >= 1L && identical(as.character(dt$Status[1]), "OK")
    }, error = function(e) FALSE)
    !ok
  }, mc.cores = n_cores))

  failed_files <- files[is_failed]
  n_ok <- sum(!is_failed)
  n_fail <- length(failed_files)
  message(sprintf("[PRUNE] Suite '%s': %d total checkpoints, %d OK (kept), %d FAILED/corrupt (%s).",
                  suite, length(files), n_ok, n_fail, if (dry_run) "dry-run" else "removing"))

  if (!dry_run && n_fail > 0L) {
    unlink(failed_files)
  }
}

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry_run" %in% args

prune_suite("synth", dry_run = dry_run)
prune_suite("real", dry_run = dry_run)
