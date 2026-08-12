# ============================================================
# main.R — AdSA-ALMD: Master Execution Script
#
# Run this script after loading your data into the ABS DataLab
# R session.  It:
#   1.  Loads the adsalmd package (installing it from source the
#       first time).
#   2.  Preprocesses the OD and adjacency matrices.
#   3.  Runs the full pipeline (SARA -> AdSA-ALMD -> refinement)
#       10 times in parallel (8 cores) via run_adsa_pipeline_parallel().
#   4.  Selects the best run and prints the Table-3 summary.
#
# The three steps this used to require running by hand — sara_iterated(),
# adsa_almd(), refine_misallocated_bgus() — are now chained for you by
# run_adsa_pipeline() (single run) / run_adsa_pipeline_parallel() (N runs).
# ============================================================
# Set this to your local clone of the adsalmd project before running.
# setwd("/path/to/AdSA-ALMD")

# ---- 0.  Install (first time only) and load the package -----
if (!requireNamespace("adsalmd", quietly = TRUE)) {
  install.packages(".", repos = NULL, type = "source")
}
library(adsalmd)

suppressPackageStartupMessages({
  library(parallel)
  library(compiler)
})
enableJIT(3)   # byte-compile all R functions for speed

# ---- 1.  Load and validate data ----------------------------
# Expects:  OD_matrix_11          — N x N numeric matrix
#           Adjacency_matrix_QLD_11 — N x N integer/logical matrix
# These objects should already exist in the DataLab workspace.

stopifnot(
  exists("OD_matrix_11"),
  exists("Adjacency_matrix_QLD_11")
)

W   <- as.matrix(OD_matrix_11)
adj <- as.matrix(Adjacency_matrix_QLD_11)

params <- AdSA_params
N      <- params$N

stopifnot(nrow(W) == N, ncol(W) == N,
          nrow(adj) == N, ncol(adj) == N)

# Precompute row / column sums ONCE — reused in every metric call
row_W <- rowSums(W)
col_W <- colSums(W)

message(sprintf("Data loaded: %d BGUs, %d total commuters",
                N, as.integer(sum(W))))

# ---- 2.  Run the full pipeline in parallel (10 runs x 8 cores) ----
message("\n--- Running AdSA-ALMD pipeline (10 independent runs) ---")

pipeline_out <- run_adsa_pipeline_parallel(
  n_runs      = params$n_runs,
  W           = W,
  adj         = adj,
  row_W       = row_W,
  col_W       = col_W,
  SC_min      = params$SC_min,
  Pop_min     = params$Pop_min,
  Pop_max     = params$Pop_max,
  Pop_tar     = params$Pop_tar,
  r           = params$penalty_exp,
  L           = params$L,
  l           = params$l,
  alpha       = params$alpha,
  beta        = params$beta,
  eps         = params$eps,
  P0          = params$P0,
  T0_samples  = params$T0_samples,
  n_cores     = params$n_cores
)

all_runs   <- pipeline_out$all_runs
run_summary <- pipeline_out$summary
best_run   <- pipeline_out$best_run

message("All runs complete.")

# ---- 3.  Report ---------------------------------------------
message(sprintf(
  "\nBest run: #%d | Base fitness = %.3f | GCI = %.3f | m = %d",
  run_summary$best_run,
  best_run$base_fitness,
  best_run$gci,
  best_run$n_clusters
))
message(sprintf(
  "Across %d runs: Mean = %.3f | SD = %.3f",
  params$n_runs,
  run_summary$mean_fitness,
  run_summary$sd_fitness
))

# ---- 4.  Misallocation count -------------------------------
message("\nCounting misallocated BGUs (this may take a moment)...")
n_misalloc <- count_misallocated(best_run$best_sol, W, adj, row_W, col_W)
message(sprintf("Misallocated BGUs: %d / %d", n_misalloc, N))

# ---- 5.  Results table (Table 3 format) --------------------
# Add other method results here for comparison, e.g.:
# results <- list(TTWA = ttwa_result, GEA = gea_result,
#                 MSA  = msa_result,  AdSA_ALMD = best_run)
results <- list(AdSA_ALMD = best_run)

cat("\n=== Results Table (Table 3) ===\n")
print(build_results_table(results, W, adj, row_W, col_W))

# ---- 6.  Save outputs --------------------------------------
saveRDS(best_run,  "outputs/best_run_AdSA_ALMD.rds")
saveRDS(all_runs,  "outputs/all_runs_AdSA_ALMD.rds")
write.csv(best_run$lma_metrics,
          "outputs/lma_metrics_best.csv", row.names = FALSE)
message("\nOutputs saved to outputs/")
