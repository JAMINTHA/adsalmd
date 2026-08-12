# ============================================================
# FILE: R/08_diagnostics.R
# Misallocation analysis and results summary (Table 3 in paper).
# ============================================================

#' Check whether a single BGU is cohesion-misallocated
#'
#' A BGU g is misallocated if there exists a neighbouring LMA such
#' that reassigning g to it would increase CI(g) without causing a
#' substantial loss in self-containment of either the origin or
#' destination market (paper: "without substantial loss in self-
#' containment").
#'
#' @param g       Integer, BGU index.
#' @param sol     Integer vector.
#' @param W       Numeric matrix.
#' @param adj     Integer/logical matrix.
#' @param row_W   Numeric vector.
#' @param col_W   Numeric vector.
#' @param sc_tol  Numeric, maximum tolerated SC reduction (default 0.05).
#' @return Logical.
is_bgu_misallocated <- function(g, sol, W, adj, row_W, col_W, sc_tol = 0.05) {
  current_cl  <- sol[g]
  current_mbr <- setdiff(cluster_members(current_cl, sol), g)
  current_ci  <- compute_ci(g, current_mbr, W, row_W, col_W)
  current_sc  <- compute_sc(cluster_members(current_cl, sol),
                             W, row_W, col_W)$sc

  adj_cls <- setdiff(unique(sol[which(adj[g, ] > 0L)]), current_cl)
  if (length(adj_cls) == 0L) return(FALSE)

  for (cl in adj_cls) {
    nbr_mbrs <- cluster_members(cl, sol)
    new_ci   <- compute_ci(g, nbr_mbrs, W, row_W, col_W)

    if (new_ci > current_ci) {
      # Would SC of the destination market drop substantially?
      candidate     <- sol
      candidate[g]  <- cl
      new_sc_dest   <- compute_sc(cluster_members(cl, candidate),
                                   W, row_W, col_W)$sc
      orig_sc_dest  <- compute_sc(nbr_mbrs, W, row_W, col_W)$sc

      if ((orig_sc_dest - new_sc_dest) <= sc_tol) return(TRUE)
    }
  }
  FALSE
}

#' Count cohesion-misallocated BGUs in a solution
#'
#' @param sol    Integer vector.
#' @param W      Numeric matrix.
#' @param adj    Integer/logical matrix.
#' @param row_W  Numeric vector.
#' @param col_W  Numeric vector.
#' @param sc_tol Numeric, tolerated SC reduction.
#' @return Integer count.
count_misallocated <- function(sol, W, adj, row_W, col_W, sc_tol = 0.05) {
  sum(vapply(seq_along(sol), is_bgu_misallocated,
             logical(1L),
             sol    = sol,
             W      = W,
             adj    = adj,
             row_W  = row_W,
             col_W  = col_W,
             sc_tol = sc_tol))
}

#' Aggregate summary statistics across multiple runs
#'
#' @param runs_list  List of result objects from adsa_almd().
#' @return data.frame with best, mean and SD of base_fitness.
summarise_runs <- function(runs_list) {
  bf <- sapply(runs_list, `[[`, "base_fitness")
  data.frame(
    best_fitness = max(bf),
    mean_fitness = mean(bf),
    sd_fitness   = sd(bf),
    best_run     = which.max(bf)
  )
}

#' Comparison table matching Table 3 in the paper
#'
#' @param results_list Named list: each element is a result from adsa_almd()
#'                     (or an equivalent structure from TTWA/GEA/MSA).
#' @param W       Numeric matrix.
#' @param adj     Integer/logical matrix.
#' @param row_W   Numeric vector.
#' @param col_W   Numeric vector.
#' @param sc_tol  Numeric, SC tolerance for misallocation.
#' @return data.frame.
build_results_table <- function(results_list, W, adj, row_W, col_W,
                                 sc_tol = 0.05) {
  do.call(rbind, lapply(names(results_list), function(method) {
    res <- results_list[[method]]
    sol <- res$best_sol
    lma <- res$lma_metrics
    data.frame(
      Method           = method,
      Base_Fitness     = round(res$base_fitness, 3),
      GCI              = round(res$gci,          3),
      SCSS             = round(mean(lma$scss),   3),
      SCDS             = round(mean(lma$scds),   3),
      m                = res$n_clusters,
      Misallocated_BGU = count_misallocated(sol, W, adj, row_W, col_W, sc_tol),
      stringsAsFactors = FALSE
    )
  }))
}
