# ============================================================
# examples/generate_random_example.R
#
# Alternative to generate_toy_examples.R's rectangular grid-block
# true partitions. That script's make_grid_true_sol() builds clusters
# as a row-band x col-band Cartesian product, so cluster sizes are
# always row_band_size * col_band_size -- you can't land on an
# arbitrary set of sizes like 2,3,4,5 with it (see conversation).
#
# This script instead grows the true clusters randomly over the
# adjacency graph (multi-source flood-fill / randomized region
# growing): pick n random seed BGUs, then repeatedly hand each
# not-yet-assigned BGU adjacent to the growing region to whichever
# cluster it touches. Every cluster is still guaranteed spatially
# contiguous by construction (each one grows outward from its seed),
# but shapes/sizes come out organic rather than rectangular.
#
# Standalone -- does not read or modify generate_toy_examples.R or
# any of its output files.
#
# Run from the project root:
#   Rscript examples/generate_random_example.R
# ============================================================

library(adsalmd)

out_dir <- "examples"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ---- grid adjacency (same construction as generate_toy_examples.R) ----

make_grid_adj <- function(nrow, ncol) {
  N <- nrow * ncol
  id <- function(r, c) (r - 1L) * ncol + c
  adj <- matrix(0L, N, N)
  for (r in seq_len(nrow)) {
    for (c in seq_len(ncol)) {
      i <- id(r, c)
      if (c < ncol) { j <- id(r, c + 1L); adj[i, j] <- 1L; adj[j, i] <- 1L }
      if (r < nrow) { j <- id(r + 1L, c); adj[i, j] <- 1L; adj[j, i] <- 1L }
    }
  }
  adj
}

# ---- random contiguous true partition (multi-source flood-fill) ----
#
# Picks k random seed BGUs, then repeatedly assigns a random
# unassigned BGU adjacent to the currently-assigned region to
# whichever neighbouring cluster it touches. Every cluster grows
# outward from its seed, so it is connected (contiguous) by
# construction -- sizes/shapes come out random, not dictated.
make_random_true_sol <- function(adj, k, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N <- nrow(adj)
  sol <- rep(0L, N)

  # sample(x, 1) samples from 1:x when x is a single number instead of
  # picking from the length-1 vector x -- guard against that here.
  pick_one <- function(x) if (length(x) == 1L) x else sample(x, 1L)

  seeds <- sample(seq_len(N), k)
  sol[seeds] <- seq_len(k)

  unassigned <- setdiff(seq_len(N), seeds)
  while (length(unassigned) > 0L) {
    assigned <- which(sol > 0L)
    frontier <- intersect(unassigned,
                           which(colSums(adj[assigned, , drop = FALSE]) > 0L))
    if (length(frontier) == 0L) {
      # disconnected leftover component -- drop it into a random cluster
      pick <- pick_one(unassigned)
      sol[pick] <- pick_one(seq_len(k))
      unassigned <- setdiff(unassigned, pick)
      next
    }
    pick <- pick_one(frontier)
    nbrs <- which(adj[pick, ] > 0L | adj[, pick] > 0L)
    nbr_assigned <- intersect(nbrs, assigned)
    sol[pick] <- sol[pick_one(nbr_assigned)]
    unassigned <- setdiff(unassigned, pick)
  }
  sol
}

# ---- synthetic W (same construction as generate_toy_examples.R) ----

make_synthetic_W <- function(true_sol, adj,
                              within_range = c(40L, 90L),
                              border_range = c(2L, 8L),
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N <- length(true_sol)
  same_cl <- outer(true_sol, true_sol, `==`)
  diag(same_cl) <- FALSE

  W <- matrix(0, N, N)
  n_same <- sum(same_cl)
  W[same_cl] <- sample(within_range[1]:within_range[2], n_same, replace = TRUE)

  border_mask <- (adj > 0) & !same_cl
  n_border <- sum(border_mask)
  if (n_border > 0) {
    W[border_mask] <- sample(border_range[1]:border_range[2], n_border, replace = TRUE)
  }
  W
}

# ---- build + save one random-partition example ----

build_random_example <- function(name, nrow, ncol, k, seed) {
  adj      <- make_grid_adj(nrow, ncol)
  true_sol <- make_random_true_sol(adj, k, seed = seed)
  stopifnot(all_clusters_contiguous(true_sol, adj))

  W        <- make_synthetic_W(true_sol, adj, seed = seed)
  row_W    <- rowSums(W)
  col_W    <- colSums(W)

  thresholds <- suggest_thresholds(true_sol, W, row_W, col_W, margin = 0.10)

  example <- list(
    name        = name,
    N           = nrow * ncol,
    grid        = c(nrow = nrow, ncol = ncol),
    W           = W,
    adj         = adj,
    row_W       = row_W,
    col_W       = col_W,
    true_sol    = true_sol,
    n_true_clusters = n_clusters(true_sol),
    thresholds  = thresholds
  )

  saveRDS(example, file.path(out_dir, paste0("toy_", name, ".rds")))
  write.csv(adj, file.path(out_dir, paste0("adj_", name, ".csv")), row.names = FALSE)
  write.csv(W, file.path(out_dir, paste0("W_", name, ".csv")), row.names = FALSE)
  write.csv(
    data.frame(bgu_id = seq_len(nrow * ncol), true_cluster = true_sol),
    file.path(out_dir, paste0("true_sol_", name, ".csv")),
    row.names = FALSE
  )

  cat(sprintf(
    "\n=== %s (N = %d, grid %dx%d, %d true clusters -- sizes: %s) ===\n",
    name, nrow * ncol, nrow, ncol, n_clusters(true_sol),
    paste(as.integer(table(true_sol)), collapse = ",")
  ))
  print(thresholds$lma_df)
  cat(sprintf(
    "Suggested: SC_min = %.3f | Pop_min = %.0f | Pop_max = %.0f | Pop_tar = %.0f\n",
    thresholds$SC_min, thresholds$Pop_min, thresholds$Pop_max, thresholds$Pop_tar
  ))

  example
}

# ---- N = 14: 2x7 grid, 4 randomly-grown contiguous clusters ----
ex14_random <- build_random_example("N14_random", nrow = 2L, ncol = 7L,
                                     k = 4L, seed = 42L)

cat("\nSaved to examples/: toy_N14_random.rds, adj_N14_random.csv,\n")
cat("                     W_N14_random.csv, true_sol_N14_random.csv\n")

# ---- verify: does run_adsa_pipeline() recover it? ----

cat("\n=== Verifying recovery on the random-partition example ===\n")
result <- run_adsa_pipeline(
  W = ex14_random$W, adj = ex14_random$adj,
  row_W = ex14_random$row_W, col_W = ex14_random$col_W,
  SC_min  = ex14_random$thresholds$SC_min,
  Pop_min = ex14_random$thresholds$Pop_min,
  Pop_max = ex14_random$thresholds$Pop_max,
  Pop_tar = ex14_random$thresholds$Pop_tar,
  r = 3L, L = 100L, l = 20L, T0_samples = 50L,
  seed = 42L, verbose = FALSE
)

check <- partitions_match(ex14_random$true_sol, result$best_sol)
cat(sprintf("Recovered clusters: %d (true: %d)\n",
            result$n_clusters, ex14_random$n_true_clusters))
cat(sprintf("Misallocated BGUs: %d / %d\n",
            count_misallocated(result$best_sol, ex14_random$W, ex14_random$adj,
                                ex14_random$row_W, ex14_random$col_W),
            ex14_random$N))
cat(sprintf("EXACT MATCH TO TRUE PARTITION: %s\n", check$match))
if (!check$match) {
  cat("Contingency table (true rows x predicted cols):\n")
  print(check$table)
}
