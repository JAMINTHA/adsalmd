# ============================================================
# examples/generate_toy_examples.R
#
# Builds three synthetic labour-market datasets with a KNOWN true
# partition baked in, for N = 14, 100, 520 BGUs. Because the true
# clusters are known by construction, you can run the pipeline and
# check it recovers them exactly (see verify_toy_examples.R) --
# this is what you'd hand to ABS as a "does the method work" sanity
# check, independent of any real confidential data.
#
# Construction:
#   - BGUs are laid out on a rectangular grid (row-major IDs).
#   - adj = rook contiguity (up/down/left/right neighbours on the grid).
#   - The true partition is a set of contiguous rectangular blocks of
#     the grid (row-bands x col-bands), so every true cluster is
#     guaranteed spatially contiguous by construction.
#   - W: strong random commuting flow between BGUs in the SAME true
#     cluster; weak "border leakage" flow only between adjacent BGUs
#     in DIFFERENT clusters (simulates realistic minor cross-border
#     commuting); zero flow between non-adjacent BGUs in different
#     clusters (no long-range noise).
#
# Run from the project root:
#   Rscript examples/generate_toy_examples.R
# ============================================================

library(adsalmd)

out_dir <- "examples"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ---- grid / adjacency helpers --------------------------------

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

# row_bands / col_bands: integer vectors of band sizes (must sum to
# nrow / ncol respectively). Returns a length-N true partition vector,
# row-major to match make_grid_adj()'s BGU numbering.
make_grid_true_sol <- function(nrow, ncol, row_bands, col_bands) {
  stopifnot(sum(row_bands) == nrow, sum(col_bands) == ncol)
  row_band_id <- rep(seq_along(row_bands), row_bands)
  col_band_id <- rep(seq_along(col_bands), col_bands)
  n_col_bands <- length(col_bands)
  cl <- matrix(0L, nrow, ncol)
  for (r in seq_len(nrow)) {
    for (c in seq_len(ncol)) {
      cl[r, c] <- (row_band_id[r] - 1L) * n_col_bands + col_band_id[c]
    }
  }
  as.vector(t(cl))
}

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

build_example <- function(name, nrow, ncol, row_bands, col_bands, seed) {
  adj      <- make_grid_adj(nrow, ncol)
  true_sol <- make_grid_true_sol(nrow, ncol, row_bands, col_bands)
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
  write.csv(adj, file.path(out_dir, paste0("adj_", name, ".csv")),
            row.names = FALSE)
  write.csv(W, file.path(out_dir, paste0("W_", name, ".csv")),
            row.names = FALSE)
  write.csv(
    data.frame(bgu_id = seq_len(nrow * ncol), true_cluster = true_sol),
    file.path(out_dir, paste0("true_sol_", name, ".csv")),
    row.names = FALSE
  )

  cat(sprintf("\n=== %s (N = %d, grid %dx%d, %d true clusters) ===\n",
              name, nrow * ncol, nrow, ncol, n_clusters(true_sol)))
  print(thresholds$lma_df)
  cat(sprintf(
    "Suggested: SC_min = %.3f | Pop_min = %.0f | Pop_max = %.0f | Pop_tar = %.0f\n",
    thresholds$SC_min, thresholds$Pop_min, thresholds$Pop_max, thresholds$Pop_tar
  ))

  example
}

# ---- N = 14: 2 rows x 7 cols, 2 clusters (cols 1-3 | cols 4-7) ----
ex14 <- build_example("N14", nrow = 2L, ncol = 7L,
                       row_bands = c(2L), col_bands = c(3L, 4L),
                       seed = 101L)

# ---- N = 100: 10x10 grid, 4 clusters (2x2 blocks of 5x5) ----
ex100 <- build_example("N100", nrow = 10L, ncol = 10L,
                        row_bands = c(5L, 5L), col_bands = c(5L, 5L),
                        seed = 102L)

# ---- N = 520: 20x26 grid, 10 clusters (5 row-bands x 2 col-bands of 4x13) ----
ex520 <- build_example("N520", nrow = 20L, ncol = 26L,
                        row_bands = c(4L, 4L, 4L, 4L, 4L), col_bands = c(13L, 13L),
                        seed = 103L)

cat("\nSaved to examples/: toy_N14.rds, toy_N100.rds, toy_N520.rds\n")
cat("                     adj_N14.csv, adj_N100.csv, adj_N520.csv\n")
cat("                     W_N14.csv, W_N100.csv, W_N520.csv\n")
cat("                     true_sol_N14.csv, true_sol_N100.csv, true_sol_N520.csv\n")

