# ============================================================
# FILE: R/04_contiguity.R
# Contiguity checking (BFS) and validity score.
# Replaces the original checkClusValidity() / is_cluster_contiguous().
# ============================================================

#' Check if a single cluster is spatially contiguous using BFS
#'
#' A cluster is contiguous if all its member BGUs are reachable
#' from any one member via adjacency links that stay within the cluster.
#' BFS is O(|members|) and avoids the fragile iteration limit in the
#' original code.
#'
#' @param members Integer vector, BGU indices in the cluster.
#' @param adj     Integer/logical matrix (N x N), adjacency matrix.
#' @return Logical.
is_cluster_contiguous <- function(members, adj) {
  if (length(members) <= 1L) return(TRUE)

  visited  <- members[1L]
  frontier <- members[1L]

  while (length(frontier) > 0L) {
    # All BGUs adjacent to the current frontier AND inside the cluster
    candidates <- which(colSums(adj[frontier, , drop = FALSE]) > 0L)
    new_nodes  <- setdiff(intersect(candidates, members), visited)
    visited    <- c(visited, new_nodes)
    frontier   <- new_nodes
  }

  length(visited) == length(members)
}

#' Check that EVERY cluster in a solution is contiguous
#'
#' Returns FALSE as soon as any cluster fails the BFS test.
#'
#' @param sol Integer vector of assignments.
#' @param adj Integer/logical matrix.
#' @return Logical.
all_clusters_contiguous <- function(sol, adj) {
  clusters <- unique_clusters(sol)
  for (cl in clusters) {
    if (!is_cluster_contiguous(cluster_members(cl, sol), adj)) {
      return(FALSE)
    }
  }
  TRUE
}

#' Validity score for a single market   [Equation 1]
#'
#' VS(Mi) = min(SC(Mi)/SC_min, 1) * min(Pop(Mi)/Pop_min, 1) * C(Mi)
#'
#' @param members Integer vector.
#' @param W       Numeric matrix.
#' @param row_W   Numeric vector.
#' @param col_W   Numeric vector.
#' @param adj     Integer/logical matrix.
#' @param SC_min  Numeric.
#' @param Pop_min Numeric.
#' @return Numeric in [0, 1].
validity_score <- function(members, W, row_W, col_W, adj, SC_min, Pop_min, Pop_max = Inf) {
  sc_val  <- compute_sc(members, W, row_W, col_W)$sc
  pop_val <- compute_population(members, row_W)
  contig  <- as.numeric(is_cluster_contiguous(members, adj))
  
  pop_score <- min(pop_val / Pop_min, 1)
  sc_score <- min(sc_val / SC_min, 1)
  sc_score * pop_score * contig
}
