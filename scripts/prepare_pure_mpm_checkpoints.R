#!/usr/bin/env Rscript
# ==============================================================================
# scripts/prepare_pure_mpm_checkpoints.R
# Removes cached task_*.rds checkpoints belonging to the methods targeted for
# Pure-MPM & audited-baseline re-execution, guarded by a one-time marker so
# interrupted reruns resume cleanly via scripts/run_remote_server.sh without
# re-deleting newly completed checkpoints.
#
# Scopes:
#   audited_all   : All 13 methods (Pure-MPM FBMS + IWSM + fixed Bag_STLSQ + deg-3 RatPoly + SR_Pool)
#   fbms6_sindy5  : All 6 FBMS methods + all 5 SINDy baselines (92,400 synth tasks)
#   fbms5_sindy5  : 5 0.025-affected FBMS methods + 5 SINDy baselines (84,000 synth tasks; skips IWSM)
#   fbms6         : All 6 FBMS methods (50,400 synth tasks)
#   fbms5         : Only the 5 methods that had 0.025 (42,000 synth tasks)
#   all           : Alias for audited_all
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0L && length(args) > idx) args[idx + 1L] else default
}

suite <- get_arg("--suite", "both")
scope <- get_arg("--scope", "audited_all")
dry_run <- "--dry_run" %in% args
force_clear <- "--force" %in% args

scope_methods <- switch(
  scope,
  "fbms5" = c(
    "BGNLM_Full_Tuned",
    "BGNLM_Hybrid_Poly",
    "BGLM_Oracle_Basis",
    "BGLM_Poly_Basis",
    "BGLM_Overcomplete"
  ),
  "fbms6" = c(
    "BGNLM_Full_Tuned",
    "BGNLM_Hybrid_Poly",
    "BGLM_Oracle_Basis",
    "BGLM_Poly_Basis",
    "BGLM_Overcomplete",
    "BGNLM_IWSM_Untuned"
  ),
  "fbms5_sindy5" = c(
    "BGNLM_Full_Tuned",
    "BGNLM_Hybrid_Poly",
    "BGLM_Oracle_Basis",
    "BGLM_Poly_Basis",
    "BGLM_Overcomplete",
    "SINDy_STLSQ_Oracle",
    "SINDy_STLSQ_Poly",
    "UQ_SINDy_Oracle",
    "UQ_SINDy_Poly",
    "SINDy_PI_Rational"
  ),
  "fbms6_sindy5" = c(
    "BGNLM_Full_Tuned",
    "BGNLM_Hybrid_Poly",
    "BGLM_Oracle_Basis",
    "BGLM_Poly_Basis",
    "BGLM_Overcomplete",
    "BGNLM_IWSM_Untuned",
    "SINDy_STLSQ_Oracle",
    "SINDy_STLSQ_Poly",
    "UQ_SINDy_Oracle",
    "UQ_SINDy_Poly",
    "SINDy_PI_Rational"
  ),
  "audited_all" = NULL,
  "all" = NULL,
  stop("Unknown --scope: ", scope, " (expected audited_all, fbms6_sindy5, fbms5_sindy5, fbms6, fbms5, all)")
)

clear_suite <- function(s) {
  manifest_path <- file.path("config", paste0("manifest_", s, ".csv"))
  ckpt_dir <- file.path("results", paste0("checkpoints_", s))
  dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)
  marker <- file.path(ckpt_dir, sprintf(".pure_mpm_cleared_%s", scope))

  if (file.exists(marker) && !force_clear) {
    message(sprintf("[SKIP] Suite '%s' already pruned for scope '%s' (marker: %s). Preserving resumed checkpoints.", s, scope, marker))
    return(invisible(NULL))
  }

  man <- fread(manifest_path)
  if (is.null(scope_methods)) {
    target_ids <- man$task_id
  } else {
    target_ids <- man[method_id %in% scope_methods, task_id]
  }
  target_files <- file.path(ckpt_dir, sprintf("task_%06d.rds", target_ids))
  existing <- target_files[file.exists(target_files)]

  message(sprintf(
    "[PREPARE PURE-MPM] Suite='%s' scope='%s': %d manifest tasks targeted, %d existing checkpoints on disk (%s).",
    s, scope, length(target_ids), length(existing), if (dry_run) "DRY-RUN" else "REMOVING"
  ))

  if (!dry_run) {
    if (length(existing) > 0L) {
      unlink(existing)
    }
    writeLines(c(
      paste0("timestamp=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      paste0("scope=", scope),
      paste0("removed_count=", length(existing))
    ), marker)
  }
}

if (suite %in% c("synth", "both")) clear_suite("synth")
if (suite %in% c("real", "both")) clear_suite("real")
