# adsalmd

Delineates local labour market areas (LMAs) from an origin-destination
commuting matrix and a spatial adjacency matrix, using:

1. **SARA** (Spatially Adaptive Random Aggregation) to build a valid initial
   partition (`sara_iterated()`).
2. **AdSA-ALMD** simulated annealing to improve it (`adsa_almd()`).
3. A deterministic misallocation-repair pass (`refine_misallocated_bgus()`).

`run_adsa_pipeline()` and `run_adsa_pipeline_parallel()` chain all three
steps so you don't have to call them by hand.

## Install

From the package directory (works fully offline — this is a local source
install, not a CRAN download):

```r
install.packages(".", repos = NULL, type = "source")
library(adsalmd)
```

`main.R` does this automatically if `adsalmd` isn't already installed.

## Inputs

Two `N x N` matrices, both indexed the same way (BGU `i` = row/col `i` in
both):

- `W` — the OD (origin-destination) matrix; `W[i, j]` = commuters from `i`
  to `j`.
- `adj` — adjacency matrix; `adj[i, j] != 0` if BGUs `i` and `j` are
  spatially adjacent.

Precompute once and reuse everywhere (every function takes these as
arguments rather than recomputing them):

```r
row_W <- rowSums(W)
col_W <- colSums(W)
```

Thresholds (`SC_min`, `Pop_min`, `Pop_max`, `Pop_tar`) and algorithm
parameters (`r`, `L`, `l`, `alpha`, `beta`, `eps`, `P0`, `T0_samples`) live
in `AdSA_params` (see `R/00_config.R`) — copy and override rather than
editing in place:

```r
params <- adsalmd::AdSA_params
params$SC_min <- 0.6
```

## One-shot run

```r
result <- run_adsa_pipeline(
  W = W, adj = adj, row_W = row_W, col_W = col_W,
  SC_min = params$SC_min, Pop_min = params$Pop_min,
  Pop_max = params$Pop_max, Pop_tar = params$Pop_tar,
  r = params$penalty_exp,
  L = params$L, l = params$l,
  seed = 42L
)

result$best_sol      # integer vector, length N — the final partition
result$base_fitness  # m * GCI, for cross-method comparison
result$gci
result$n_clusters
result$lma_metrics   # per-market data.frame: cluster, pop, scss, scds, sc, n_bgus
```

## Parallel run (multiple independent attempts)

Runs the same pipeline `n_runs` times on a local cluster (different random
seeds) and returns every run plus the best one by `base_fitness`:

```r
out <- run_adsa_pipeline_parallel(
  n_runs = 10L,
  W = W, adj = adj, row_W = row_W, col_W = col_W,
  SC_min = params$SC_min, Pop_min = params$Pop_min,
  Pop_max = params$Pop_max, Pop_tar = params$Pop_tar,
  r = params$penalty_exp,
  L = params$L, l = params$l,
  n_cores = params$n_cores
)

out$best_run     # the single best result (same structure as run_adsa_pipeline())
out$summary      # data.frame: best_fitness, mean_fitness, sd_fitness, best_run (index)
out$all_runs     # list of all n_runs results, for inspecting spread
```

Each worker loads the package with `library(adsalmd)`, so the package must
be installed (not just source()'d) before calling this.

## Inspecting / validating a solution

Given any `sol` (an integer vector of length N, e.g. `result$best_sol`):

| Function | What it tells you |
|---|---|
| `build_lma_df(sol, W, row_W, col_W)` | Per-market table: population, SCSS, SCDS, SC, BGU count |
| `is_cluster_contiguous(members, adj)` / `all_clusters_contiguous(sol, adj)` | Whether a market (or every market) is spatially contiguous |
| `compute_gci(sol, W, row_W, col_W)` / `compute_base_fitness(sol, W, row_W, col_W)` | Global cohesion index / `m * GCI` |
| `compute_sc()`, `compute_scss()`, `compute_scds()` | Self-containment components for a set of members |
| `validity_score(members, W, row_W, col_W, adj, SC_min, Pop_min)` | Eq. 1 validity score for one market |
| `is_bgu_misallocated(g, sol, W, adj, row_W, col_W)` / `count_misallocated(sol, W, adj, row_W, col_W)` | Whether a single BGU (or how many BGUs) would fit better in a neighbouring market |
| `summarise_runs(all_runs)` | Best/mean/SD of `base_fitness` across a list of runs |
| `build_results_table(list(Method = result), W, adj, row_W, col_W)` | Paper-style Table 3 comparison row(s) |

Every exported function has a help page — `?run_adsa_pipeline`,
`?build_lma_df`, etc. (roxygen2 isn't installable on this R version, so
these `.Rd` files under `man/` are hand-maintained rather than
auto-generated; if you get access to a newer R + roxygen2 later, running
`roxygen2::roxygenise()` will regenerate them from the `#'` comments
already in the `R/` source files).

## `main.R`

`main.R` is the DataLab entry-point script: loads `OD_matrix_11` /
`Adjacency_matrix_QLD_11` from the workspace, installs/loads `adsalmd`,
calls `run_adsa_pipeline_parallel()` with `AdSA_params`, prints the
Table-3 summary, and saves `outputs/best_run_AdSA_ALMD.rds`,
`outputs/all_runs_AdSA_ALMD.rds`, `outputs/lma_metrics_best.csv`.

## Known rough edges

- `sara_iterated()` prints progress unconditionally (its `verbose` argument
  isn't actually wired to anything internally yet).
- `operator_1`, `operator_6`, `operator_8`, `operator_10`, `operator_11`
  (and the `apply_random_operator()` dispatcher) read a global `adj`
  matrix and a global `params` list ($SC_min, $Pop_min, $Pop_max) instead
  of taking them as arguments. `run_adsa_pipeline()` /
  `run_adsa_pipeline_parallel()` set these up for you (and restore
  whatever was there before), so you don't need to worry about it unless
  you're calling `adsa_almd()` or the operators directly outside the
  pipeline — in that case, set `adj` and `params` as globals yourself
  first.
- `tests/test_functions.R` is a standalone sanity-check script (not part
  of the installed package — see `.Rbuildignore`); it still works the same
  way it always has via `source("tests/test_functions.R")`.
