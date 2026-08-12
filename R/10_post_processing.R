#' Local search refinement: reassign misallocated BGUs to a better-fitting
#' adjacent cluster.
#'
#' Repeatedly scans every BGU. If a BGU is misallocated (per
#' is_bgu_misallocated), it is reassigned to whichever adjacent cluster gives
#' the highest compute_ci score, as long as the resulting partition keeps
#' every cluster contiguous. Passes repeat until a full pass makes no further
#' improvement, or max_iter passes have been made.
#'
#' @param best_sol Integer vector, starting partition to improve.
#' @param W        Numeric matrix (N x N).
#' @param adj      Integer/logical matrix (N x N).
#' @param row_W    Numeric vector, rowSums(W).
#' @param col_W    Numeric vector, colSums(W).
#' @param max_iter Integer, maximum number of full passes over all BGUs (default 100).
#' @return         Integer vector, the improved (relabelled) partition.
refine_misallocated_bgus <- function(best_sol, W, adj, row_W, col_W, max_iter = 100L) {
  
  improved <- TRUE
  iter <- 0L
  
  while (improved && iter < max_iter) {
    iter <- iter + 1L
    improved <- FALSE
    
    for (g in seq_along(best_sol)) {
      if (is_bgu_misallocated(g, best_sol, W, adj, row_W, col_W)) {
        
        cur_cl <- best_sol[g]
        adj_cls <- setdiff(
          unique(best_sol[which(adj[g, ] > 0L)]), cur_cl
        )
        
        # nothing to move to
        if (length(adj_cls) == 0L) next
        
        scores <- vapply(adj_cls, function(cl)
          compute_ci(g, cluster_members(cl, best_sol), W, row_W, col_W),
          numeric(1L)
        )
        
        best_cl <- adj_cls[which.max(scores)]
        cand <- best_sol
        cand[g] <- best_cl
        cand <- getOrderedVec(cand)
        
        if (all_clusters_contiguous(cand, adj)) {
          best_sol <- cand
          improved <- TRUE
        }
      }
    }
  }
  
  best_sol
}