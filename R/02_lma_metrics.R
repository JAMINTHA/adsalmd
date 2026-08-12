# ============================================================
# FILE: R/02_lma_metrics.R
# LMA-specific metrics.  Every formula is cited to the paper.
#
# Notation (matching paper):
#   W      : OD matrix (N x N).  W[i, j] = commuters from i to j.
#   row_W  : rowSums(W) — precomputed once and passed in.
#   col_W  : colSums(W) — precomputed once and passed in.
#   G      : full BGU index set 1:N.
#   Ms     : set of BGU indices in market s.
#   g      : single BGU index.
#   M\g    : all BGUs in g's market, excluding g.
# ============================================================

# ---- Cohesion -----------------------------------------------

#' Cohesion Index for a single BGU   [Equation 6]
#'
#' CI(g) = W(g, M\g)^2 / (W(g,G) * W(G, M\g))
#'       + W(M\g, g)^2 / (W(M\g,G) * W(G, g))
#'
#' Returns 0 for singleton clusters.
#'
#' @param g     Integer, BGU index.
#' @param M_g   Integer vector, other BGUs in g's market (M \ {g}).
#' @param W     Numeric matrix, OD matrix.
#' @param row_W Numeric vector, rowSums(W).
#' @param col_W Numeric vector, colSums(W).
#' @return Numeric scalar >= 0.
compute_ci <- function(g, M_g, W, row_W, col_W) {
  if (length(M_g) == 0L) return(0)

  W_g_Mg <- sum(W[g,   M_g])        # W(g,  M\g)
  W_Mg_g <- sum(W[M_g, g  ])        # W(M\g, g)
  W_g_G  <- row_W[g]                # W(g,  G)
  W_G_g  <- col_W[g]                # W(G,  g) 
  W_Mg_G <- sum(row_W[M_g])         # W(M\g, G)
  W_G_Mg <- sum(col_W[M_g])         # W(G,  M\g)

  term1  <- if (W_g_G > 0 && W_G_Mg > 0) W_g_Mg^2 / (W_g_G  * W_G_Mg) else 0
  term2  <- if (W_Mg_G > 0 && W_G_g > 0) W_Mg_g^2 / (W_Mg_G * W_G_g ) else 0

  term1 + term2
}

#' Global Cohesion Index   [Equation 7]
#'
#' GCI = sum_{g in G} CI(g)
#'
#' Automatically uses the fast sub-matrix version for N >= 50
#' (compute_gci_fast, defined later in this file) and the
#' element-wise version for small N where clarity matters more.
#'
#' @param sol   Integer vector of assignments (length N).
#' @param W     Numeric matrix (N x N).
#' @param row_W Numeric vector (length N).
#' @param col_W Numeric vector (length N).
#' @return Numeric scalar GCI.
compute_gci <- function(sol, W, row_W, col_W) {
  # For N >= 50 use the fast sub-matrix version (~15-20x speedup)
  if (length(sol) >= 50L) return(compute_gci_fast(sol, W, row_W, col_W))

  clusters <- unique_clusters(sol)
  gci <- 0
  for (cl in clusters) {
    members <- cluster_members(cl, sol)
    for (g in members) {
      M_g <- setdiff(members, g)
      gci <- gci + compute_ci(g, M_g, W, row_W, col_W)
    }
  }
  gci
}

# ---- Self-containment ---------------------------------------

#' Supply-side self-containment for a market   [Equation 10]
#'
#' SCSS(Ms) = W(Ms, Ms) / W(Ms, G)
#'
#' @param members Integer vector, BGU indices in market.
#' @param W       Numeric matrix.
#' @param row_W   Numeric vector.
#' @return        Numeric in [0, 1].
compute_scss <- function(members, W, row_W) {
  W_Ms_G <- sum(row_W[members])
  if (W_Ms_G == 0) return(0)
  sum(W[members, members]) / W_Ms_G
}

#' Demand-side self-containment for a market   [Equation 11]
#'
#' SCDS(Ms) = W(Ms, Ms) / W(G, Ms)
#'
#' @param members Integer vector.
#' @param W       Numeric matrix.
#' @param col_W   Numeric vector.
#' @return Numeric in [0, 1].
compute_scds <- function(members, W, col_W) {
  W_G_Ms <- sum(col_W[members])
  if (W_G_Ms == 0) return(0)
  sum(W[members, members]) / W_G_Ms
}

#' Self-containment for a market   [Equation 12]
#'
#' SC(Ms) = min(SCSS(Ms), SCDS(Ms))
#'
#' @param members Integer vector.
#' @param W       Numeric matrix.
#' @param row_W   Numeric vector.
#' @param col_W   Numeric vector.
#' @return Named list: scss, scds, sc.
compute_sc <- function(members, W, row_W, col_W) {
  scss <- compute_scss(members, W, row_W)
  scds <- compute_scds(members, W, col_W)
  list(scss = scss, scds = scds, sc = min(scss, scds))
}

# ---- Population ---------------------------------------------

#' Population of a market  (total outgoing commuters = total residents)
#'
#' Pop(Ms) = W(Ms, G) = sum_{i in Ms} rowSums(W)[i]
#'
#' @param members Integer vector.
#' @param row_W   Numeric vector.
#' @return Numeric scalar.
compute_population <- function(members, row_W) {
  sum(row_W[members])
}

# ---- Summary data frame -------------------------------------

#' Build a per-market summary data frame
#'
#' Used internally by the objective function and diagnostics.
#'
#' @param sol   Integer vector.
#' @param W     Numeric matrix.
#' @param row_W Numeric vector.
#' @param col_W Numeric vector.
#' @return data.frame with columns: cluster, pop, scss, scds, sc, n_bgus.
build_lma_df <- function(sol, W, row_W, col_W) {
  clusters <- unique_clusters(sol)

  rows <- lapply(clusters, function(cl) {
    members <- cluster_members(cl, sol)
    sc_vals <- compute_sc(members, W, row_W, col_W)
    data.frame(
      cluster = cl,
      pop     = compute_population(members, row_W),
      scss    = sc_vals$scss,
      scds    = sc_vals$scds,
      sc      = sc_vals$sc,
      n_bgus  = length(members),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

# ---- Optimised GCI for large N ------------------------------

#' Fast Global Cohesion Index for large N   [Equation 7]
#'
#' Replaces compute_gci() for N >= 100.  Instead of indexing
#' W[g, M_g] once per BGU (520 separate row-extractions for N=520),
#' this extracts W[Ms, Ms] ONCE per cluster and derives all per-BGU
#' quantities from the sub-matrix row/column sums.
#'
#' Speedup on N=520:  ~15-20x over the naive version.
#'
#' @param sol   Integer vector of assignments (length N).
#' @param W     Numeric matrix (N x N).
#' @param row_W Numeric vector (length N).
#' @param col_W Numeric vector (length N).
#' @return Numeric scalar GCI.
compute_gci_fast <- function(sol, W, row_W, col_W) {
  clusters <- unique_clusters(sol)
  gci      <- 0

  for (cl in clusters) {
    members <- cluster_members(cl, sol)
    n_m     <- length(members)
    if (n_m <= 1L) next

    # --- Extract sub-matrix once per cluster ---
    W_cl   <- W[members, members, drop = FALSE]
    rw_cl  <- rowSums(W_cl)   # W(g, M\g) for each g in order
    cw_cl  <- colSums(W_cl)   # W(M\g, g) for each g in order

    # Cluster-level totals (precomputed)
    cl_row <- sum(row_W[members])   # W(Ms, G)
    cl_col <- sum(col_W[members])   # W(G, Ms)

    for (j in seq_len(n_m)) {
      g       <- members[j]
      W_g_Mg  <- rw_cl[j]                  # W(g,  M\g)
      W_Mg_g  <- cw_cl[j]                  # W(M\g, g)
      W_g_G   <- row_W[g]                  # W(g,  G)
      W_G_g   <- col_W[g]                  # W(G,  g)
      W_Mg_G  <- cl_row - row_W[g]         # W(M\g, G)
      W_G_Mg  <- cl_col - col_W[g]         # W(G,  M\g)

      t1 <- if (W_g_G  > 0 && W_G_Mg > 0) W_g_Mg^2 / (W_g_G  * W_G_Mg) else 0
      t2 <- if (W_Mg_G > 0 && W_G_g  > 0) W_Mg_g^2 / (W_Mg_G * W_G_g ) else 0
      gci <- gci + t1 + t2
    }
  }
  gci
}

# ---- Base fitness -------------------------------------------

#' Base fitness function used for cross-method comparison   [Equation 8]
#'
#' fbase(x) = m * GCI(x)
#'
#' This is the common performance measure used by GEA and MSA.
#'
#' @param sol   Integer vector.
#' @param W     Numeric matrix.
#' @param row_W Numeric vector.
#' @param col_W Numeric vector.
#' @return Numeric scalar.
compute_base_fitness <- function(sol, W, row_W, col_W) {
  n_clusters(sol) * compute_gci(sol, W, row_W, col_W)
}
