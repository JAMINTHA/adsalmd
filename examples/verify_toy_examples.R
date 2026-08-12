# ============================================================
# examples/verify_toy_examples.R
#
# Runs run_adsa_pipeline() on each of the three toy datasets
# (generate_toy_examples.R must be run first) and checks that the
# recovered partition matches the known true partition exactly
# (up to arbitrary cluster relabelling), via partitions_match().
#
# L/l here are kept small since these are easy, cleanly-separable
# synthetic cases -- for real data use the full AdSA_params defaults.
# ============================================================

library(adsalmd)

verify_one <- function(path, L = 100L, l = 20L, T0_samples = 50L, seed = 42L) {
  ex <- readRDS(path)
  cat(sprintf("\n=== %s (N = %d, %d true clusters) ===\n",
              ex$name, ex$N, ex$n_true_clusters))

  result <- run_adsa_pipeline(
    W = ex$W, adj = ex$adj, row_W = ex$row_W, col_W = ex$col_W,
    SC_min  = ex$thresholds$SC_min,
    Pop_min = ex$thresholds$Pop_min,
    Pop_max = ex$thresholds$Pop_max,
    Pop_tar = ex$thresholds$Pop_tar,
    r = 3L, L = L, l = l, T0_samples = T0_samples,
    seed = seed, verbose = FALSE
  )

  check <- partitions_match(ex$true_sol, result$best_sol)

  cat(sprintf("Recovered clusters: %d (true: %d)\n",
              result$n_clusters, ex$n_true_clusters))
  cat(sprintf("Misallocated BGUs: %d / %d\n",
              count_misallocated(result$best_sol, ex$W, ex$adj, ex$row_W, ex$col_W),
              ex$N))
  cat(sprintf("EXACT MATCH TO TRUE PARTITION: %s\n", check$match))
  if (!check$match) {
    cat("Contingency table (true rows x predicted cols):\n")
    print(check$table)
  }

  list(result = result, check = check)
}

r14  <- verify_one("examples/toy_N14.rds")
r100 <- verify_one("examples/toy_N100.rds", L = 150L, l = 20L)
r520 <- verify_one("examples/toy_N520.rds", L = 200L, l = 25L, T0_samples = 100L)

cat("\n=== SUMMARY ===\n")
cat("N14  match:", r14$check$match,  "\n")
cat("N100 match:", r100$check$match, "\n")
cat("N520 match:", r520$check$match, "\n")
