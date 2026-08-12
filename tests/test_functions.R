# ============================================================
# tests/test_functions.R
# Unit tests — run this to verify correctness BEFORE using
# the algorithm on real data.
#
# No external testing framework needed; uses base R stopifnot().
# Run with:  source("tests/test_functions.R")
# ============================================================

source("R/00_config.R"); source("R/01_utils.R")
source("R/02_lma_metrics.R"); source("R/03_objective.R")
source("R/04_contiguity.R"); source("R/05_sara.R")
source("R/06_operators.R"); source("R/07_adsa_almd.R")
source("R/08_diagnostics.R")

cat("Running unit tests...\n")

# ============================================================
# Tiny synthetic dataset: 6 BGUs in 2 LMAs
#   LMA 1: BGUs 1, 2, 3  (arranged in a line: 1-2-3)
#   LMA 2: BGUs 4, 5, 6  (arranged in a line: 4-5-6)
#   BGUs 3 and 4 are adjacent (border between LMAs)
# ============================================================

N_test <- 6L

# OD matrix: strong internal flows, weak cross-LMA flows
set.seed(1L)
W_test <- matrix(0, N_test, N_test)
W_test[1, 2] <- 80; W_test[2, 1] <- 70
W_test[2, 3] <- 90; W_test[3, 2] <- 85
W_test[1, 3] <- 60; W_test[3, 1] <- 55
W_test[4, 5] <- 80; W_test[5, 4] <- 75
W_test[5, 6] <- 85; W_test[6, 5] <- 80
W_test[4, 6] <- 65; W_test[6, 4] <- 60
# Cross-LMA (weak)
W_test[3, 4] <- 10; W_test[4, 3] <- 8

# Adjacency: 1-2-3-4-5-6 in a chain
adj_test <- matrix(0L, N_test, N_test)
adj_test[1, 2] <- adj_test[2, 1] <- 1L
adj_test[2, 3] <- adj_test[3, 2] <- 1L
adj_test[3, 4] <- adj_test[4, 3] <- 1L
adj_test[4, 5] <- adj_test[5, 4] <- 1L
adj_test[5, 6] <- adj_test[6, 5] <- 1L

row_W_test <- rowSums(W_test)
col_W_test <- colSums(W_test)

sol_true <- c(1L, 1L, 1L, 2L, 2L, 2L)   # correct partition

# ---- Test 1: get_ordered_vec --------------------------------
cat("  [1] get_ordered_vec ... ")
x <- c(3L, 3L, 7L, 7L, 1L)
stopifnot(identical(get_ordered_vec(x), c(2L, 2L, 3L, 3L, 1L)))
# Already consecutive — unchanged
stopifnot(identical(get_ordered_vec(c(1L, 2L, 3L)), c(1L, 2L, 3L)))
cat("PASS\n")

# ---- Test 2: cluster_members --------------------------------
cat("  [2] cluster_members ... ")
stopifnot(identical(cluster_members(1L, sol_true), c(1L, 2L, 3L)))
stopifnot(identical(cluster_members(2L, sol_true), c(4L, 5L, 6L)))
cat("PASS\n")

# ---- Test 3: compute_scss / scds / sc ----------------------
cat("  [3] self-containment metrics ... ")
mbrs1 <- cluster_members(1L, sol_true)
sc1   <- compute_sc(mbrs1, W_test, row_W_test, col_W_test)
# Internal flow of LMA 1 = W_test[1:3, 1:3]
internal <- sum(W_test[mbrs1, mbrs1])
scss_exp <- internal / sum(row_W_test[mbrs1])
scds_exp <- internal / sum(col_W_test[mbrs1])
stopifnot(abs(sc1$scss - scss_exp) < 1e-10)
stopifnot(abs(sc1$scds - scds_exp) < 1e-10)
stopifnot(sc1$sc == min(scss_exp, scds_exp))
cat("PASS\n")

# ---- Test 4: compute_ci / compute_gci ----------------------
cat("  [4] CI and GCI ... ")
# BGU 2's cluster-mates are {1, 3}
ci_2 <- compute_ci(2L, c(1L, 3L), W_test, row_W_test, col_W_test)
stopifnot(is.finite(ci_2) && ci_2 >= 0)
gci  <- compute_gci(sol_true, W_test, row_W_test, col_W_test)
stopifnot(is.finite(gci) && gci > 0)
# GCI must equal sum of individual CI values
ci_manual <- sum(vapply(seq_len(N_test), function(g) {
  cl   <- sol_true[g]
  M_g  <- setdiff(cluster_members(cl, sol_true), g)
  compute_ci(g, M_g, W_test, row_W_test, col_W_test)
}, numeric(1L)))
stopifnot(abs(gci - ci_manual) < 1e-10)
cat("PASS\n")

# ---- Test 5: penalty functions ------------------------------
cat("  [5] penalty functions ... ")
# No penalty if SC above threshold
stopifnot(penalty_sc_single(0.8, SC_min = 0.7, r = 3L) == 0)
# Penalty when below threshold
p <- penalty_sc_single(0.5, SC_min = 0.7, r = 3L)
stopifnot(p > 0 && p <= 1)
# Exact value: 1 - (0.5/0.7)^3
stopifnot(abs(p - (1 - (0.5 / 0.7)^3)) < 1e-12)
# Same logic for population
stopifnot(penalty_pop_single(600, Pop_min = 500, r = 3L) == 0)
p2 <- penalty_pop_single(200, Pop_min = 500, r = 3L)
stopifnot(abs(p2 - (1 - (200 / 500)^3)) < 1e-12)
cat("PASS\n")

# ---- Test 6: contiguity checks ------------------------------
cat("  [6] contiguity ... ")
stopifnot(is_cluster_contiguous(c(1L, 2L, 3L), adj_test))
stopifnot(!is_cluster_contiguous(c(1L, 3L), adj_test))   # 1 and 3 not adjacent
stopifnot(all_clusters_contiguous(sol_true, adj_test))
# Non-contiguous partition
bad_sol <- c(1L, 2L, 1L, 2L, 2L, 2L)   # LMA 1 = {1,3}, not adjacent
stopifnot(!all_clusters_contiguous(bad_sol, adj_test))
cat("PASS\n")

# ---- Test 7: validity score ---------------------------------
cat("  [7] validity_score ... ")
vs <- validity_score(cluster_members(1L, sol_true),
                     W_test, row_W_test, col_W_test,
                     adj_test, SC_min = 0.3, Pop_min = 1)
stopifnot(vs >= 0 && vs <= 1)
cat("PASS\n")

# ---- Test 8: compute_objective ------------------------------
cat("  [8] compute_objective ... ")
F_val <- compute_objective(sol_true, W_test, row_W_test, col_W_test,
                            SC_min = 0.3, Pop_min = 1,
                            r = 3L, adj = adj_test)
stopifnot(is.finite(F_val))
# Non-contiguous solution must return 0
F_bad <- compute_objective(bad_sol, W_test, row_W_test, col_W_test,
                            SC_min = 0.3, Pop_min = 1,
                            r = 3L, adj = adj_test)
stopifnot(F_bad == 0)
cat("PASS\n")

# ---- Test 9: compute_base_fitness ---------------------------
cat("  [9] compute_base_fitness ... ")
bf <- compute_base_fitness(sol_true, W_test, row_W_test, col_W_test)
stopifnot(abs(bf - 2 * gci) < 1e-10)   # m=2 clusters
cat("PASS\n")

# ---- Test 10: SARA on atomic start --------------------------
cat("  [10] SARA ... ")
set.seed(99L)
sara_sol <- sara(seq_len(N_test), W_test, adj_test,
                 row_W_test, col_W_test,
                 SC_min = 0.1, Pop_min = 1)
stopifnot(length(sara_sol) == N_test)
stopifnot(all_clusters_contiguous(sara_sol, adj_test))
stopifnot(all(sara_sol >= 1L))
cat("PASS\n")

# ---- Test 11: operator smoke test ---------------------------
cat("  [11] operators (smoke test) ... ")
for (op_fn in list(operator_1, operator_2, operator_3, operator_4,
                    operator_5, operator_6, operator_7, operator_8,
                    operator_9, operator_10)) {
  set.seed(7L)
  out <- op_fn(sol_true, W_test, adj_test,
               row_W_test, col_W_test, SC_min = 0.3, Pop_min = 1)
  stopifnot(length(out) == N_test)
  stopifnot(all(out >= 1L))
  # Any returned solution must be contiguous
  stopifnot(all_clusters_contiguous(out, adj_test))
}
cat("PASS\n")

# ---- Test 12: T0 estimation ---------------------------------
cat("  [12] estimate_t0 ... ")
set.seed(42L)
T0 <- estimate_t0(sol_true, W_test, adj_test, row_W_test, col_W_test,
                   SC_min = 0.3, Pop_min = 1, r = 3L,
                   P0 = 0.8, n_samples = 20L)
stopifnot(is.finite(T0) && T0 > 0)
cat("PASS\n")

# ---- Test 13: adsa_almd mini-run ----------------------------
cat("  [13] adsa_almd (mini-run, L=10) ... ")
set.seed(1L)
res <- adsa_almd(
  init_sol   = sol_true,
  W          = W_test, adj = adj_test,
  row_W      = row_W_test, col_W = col_W_test,
  SC_min     = 0.3, Pop_min = 1, r = 3L,
  L = 10L, l = 5L, verbose = FALSE
)
stopifnot(is.list(res))
stopifnot(all(c("best_sol", "best_fitness", "gci",
                "base_fitness", "n_clusters",
                "lma_metrics", "history") %in% names(res)))
stopifnot(all_clusters_contiguous(res$best_sol, adj_test))
cat("PASS\n")

# ---- Test 14: misallocation count ---------------------------
cat("  [14] count_misallocated ... ")
n_m <- count_misallocated(sol_true, W_test, adj_test,
                           row_W_test, col_W_test)
stopifnot(is.integer(n_m) || is.numeric(n_m))
stopifnot(n_m >= 0 && n_m <= N_test)
cat("PASS\n")

cat("\nAll 14 tests PASSED.\n")
